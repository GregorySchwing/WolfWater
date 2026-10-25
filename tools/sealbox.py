"""A libsodium sealed box in Python's standard library alone.

GitHub stores a secret only as a sealed box made for its public key
(crypto_box_seal: an X25519 key made for the one message, XSalsa20 and
Poly1305), and what made it cannot open it again. The token workflow
(.github/workflows/token.yml) seals the token it mints for the key the
control panel hands it, so the token leaves the runner only in the form
GitHub's secrets API takes, and the panel stores it without ever holding
it. PyNaCl would do this in one line; the project takes no new dependency
(.claude/rules/writing.md), so the primitives are written out here from
RFC 7748 (X25519), the Salsa20 specification (XSalsa20 and HSalsa20) and
RFC 8439 (Poly1305). open_sealed is for the tests, which check both
directions against a box libsodium made.

  python3 tools/sealbox.py seal KEY_BASE64 < message   the box, in base64
"""

import base64
import hashlib
import hmac
import os
import sys

P = 2 ** 255 - 19
A24 = 121665
MASK32 = 0xFFFFFFFF
SIGMA = b"expand 32-byte k"
KEY_BYTES = 32
TAG_BYTES = 16
NONCE_BYTES = 24


class SealError(Exception):
    """A key or a box that is not what a sealed box needs."""


def x25519(scalar, point):
    """RFC 7748's Montgomery ladder: scalar times the point's u."""
    k = bytearray(scalar)
    k[0] &= 248
    k[31] &= 127
    k[31] |= 64
    k = int.from_bytes(k, "little")
    u = int.from_bytes(point, "little") & ((1 << 255) - 1)
    x2, z2, x3, z3, swap = 1, 0, u, 1, 0
    for t in reversed(range(255)):
        bit = (k >> t) & 1
        swap ^= bit
        if swap:
            x2, x3, z2, z3 = x3, x2, z3, z2
        swap = bit
        a, b = (x2 + z2) % P, (x2 - z2) % P
        aa, bb = a * a % P, b * b % P
        e = (aa - bb) % P
        c, d = (x3 + z3) % P, (x3 - z3) % P
        da, cb = d * a % P, c * b % P
        x3 = (da + cb) ** 2 % P
        z3 = u * (da - cb) ** 2 % P
        x2 = aa * bb % P
        z2 = e * (aa + A24 * e) % P
    if swap:
        x2, z2 = x3, z3
    return (x2 * pow(z2, P - 2, P) % P).to_bytes(32, "little")


def public_key(secret):
    return x25519(secret, (9).to_bytes(32, "little"))


def rotate(value, count):
    return ((value << count) | (value >> (32 - count))) & MASK32


def quarter(s, a, b, c, d):
    s[b] ^= rotate((s[a] + s[d]) & MASK32, 7)
    s[c] ^= rotate((s[b] + s[a]) & MASK32, 9)
    s[d] ^= rotate((s[c] + s[b]) & MASK32, 13)
    s[a] ^= rotate((s[d] + s[c]) & MASK32, 18)


def salsa_rounds(state):
    s = list(state)
    for _ in range(10):
        quarter(s, 0, 4, 8, 12)
        quarter(s, 5, 9, 13, 1)
        quarter(s, 10, 14, 2, 6)
        quarter(s, 15, 3, 7, 11)
        quarter(s, 0, 1, 2, 3)
        quarter(s, 5, 6, 7, 4)
        quarter(s, 10, 11, 8, 9)
        quarter(s, 15, 12, 13, 14)
    return s


def words(data):
    return [int.from_bytes(data[i:i + 4], "little") for i in range(0, len(data), 4)]


def salsa_state(key, middle):
    """Salsa20's sixteen words: the constant on the diagonal, the key's
    halves around it, sixteen bytes of nonce and counter in the middle."""
    c, k, m = words(SIGMA), words(key), words(middle)
    return [c[0], *k[0:4], c[1], *m, c[2], *k[4:8], c[3]]


def hsalsa20(key, nonce16):
    s = salsa_rounds(salsa_state(key, nonce16))
    out = [s[i] for i in (0, 5, 10, 15, 6, 7, 8, 9)]
    return b"".join(w.to_bytes(4, "little") for w in out)


def salsa20_stream(key, nonce8, length):
    blocks = []
    for counter in range((length + 63) // 64):
        state = salsa_state(key, nonce8 + counter.to_bytes(8, "little"))
        mixed = salsa_rounds(state)
        blocks.append(b"".join(((a + b) & MASK32).to_bytes(4, "little") for a, b in zip(mixed, state)))
    return b"".join(blocks)[:length]


def poly1305(key, message):
    r = int.from_bytes(key[:16], "little") & 0x0FFFFFFC0FFFFFFC0FFFFFFC0FFFFFFF
    s = int.from_bytes(key[16:32], "little")
    p, acc = (1 << 130) - 5, 0
    for i in range(0, len(message), 16):
        block = message[i:i + 16] + b"\x01"
        acc = (acc + int.from_bytes(block, "little")) * r % p
    return ((acc + s) & ((1 << 128) - 1)).to_bytes(16, "little")


def xsalsa20_poly1305(key, nonce, message, opening=False):
    """crypto_secretbox: the first 32 bytes of the stream key Poly1305, the
    rest encrypt; opening checks the tag on the ciphertext first."""
    subkey = hsalsa20(key, nonce[:16])
    stream = salsa20_stream(subkey, nonce[16:24], 32 + len(message) - (TAG_BYTES if opening else 0))
    if opening:
        tag, body = message[:TAG_BYTES], message[TAG_BYTES:]
        if not hmac.compare_digest(tag, poly1305(stream[:32], body)):
            raise SealError("the box does not open with this key")
        return bytes(a ^ b for a, b in zip(body, stream[32:]))
    body = bytes(a ^ b for a, b in zip(message, stream[32:]))
    return poly1305(stream[:32], body) + body


def box_key(secret, peer):
    """crypto_box's shared key: HSalsa20 of the X25519 product."""
    shared = x25519(secret, peer)
    if shared == bytes(32):
        raise SealError("that public key gives no shared key")
    return hsalsa20(shared, bytes(16))


def seal_nonce(ephemeral, recipient):
    return hashlib.blake2b(ephemeral + recipient, digest_size=NONCE_BYTES).digest()


def seal(message, recipient, ephemeral_secret=None):
    """crypto_box_seal: the ephemeral public key, then the box."""
    if len(recipient) != KEY_BYTES:
        raise SealError("a sealed box's key is 32 bytes")
    secret = ephemeral_secret or os.urandom(KEY_BYTES)
    ephemeral = public_key(secret)
    box = xsalsa20_poly1305(box_key(secret, recipient), seal_nonce(ephemeral, recipient), message)
    return ephemeral + box


def open_sealed(sealed, recipient_secret):
    """crypto_box_seal_open, for the tests."""
    if len(sealed) < KEY_BYTES + TAG_BYTES:
        raise SealError("too short to be a sealed box")
    ephemeral, box = sealed[:KEY_BYTES], sealed[KEY_BYTES:]
    recipient = public_key(recipient_secret)
    return xsalsa20_poly1305(box_key(recipient_secret, ephemeral), seal_nonce(ephemeral, recipient), box, opening=True)


def key_from_base64(text):
    try:
        key = base64.b64decode(text.strip(), validate=True)
    except ValueError:
        raise SealError("the key is not base64") from None
    if len(key) != KEY_BYTES:
        raise SealError("a sealed box's key is 32 bytes")
    return key


def main(argv):
    if len(argv) != 2 or argv[0] != "seal":
        print(__doc__.strip().splitlines()[-1].strip(), file=sys.stderr)
        return 2
    try:
        sealed = seal(sys.stdin.buffer.read(), key_from_base64(argv[1]))
    except SealError as error:
        print(f"sealbox: {error}", file=sys.stderr)
        return 1
    print(base64.b64encode(sealed).decode())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

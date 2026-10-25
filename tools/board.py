"""Render the board and nudge the pull requests that wait: a shim.

The board's code now lives in the claim protocol's library,
OpenSpineConsortium/symphony-protocol, vendored here as tools/sdlc_protocol
(tools/VENDORED names the version and the commit). This file stays because
.github/workflows/board.yml and every copy of the template already made
run "python3 tools/board.py ...": it hands the same arguments to
"python3 -m sdlc_protocol board", which takes exactly the arguments and
gives exactly the exit codes this tool always did.

A comment the run posts and GitHub refuses, a "main moved" or an idle
note on a pull request, is a line of the workflow's log, and the run
still exits 0, so the board is written to its branch whatever GitHub
answered; a refused note, idle or "main moved", is posted by a later run,
since both depend on the state and not on the run's event (from the
library's 0.5.0; before it a refused "main moved" waited for the default
branch to move again).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sdlc_protocol.cli import main as protocol_main  # noqa: E402


def main(argv=None):
    return protocol_main(["board", *(sys.argv[1:] if argv is None else argv)])


if __name__ == "__main__":
    sys.exit(main())

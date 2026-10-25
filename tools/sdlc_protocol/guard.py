"""Refuse a git push the board rules forbid, before it runs.

Claude Code runs this as a PreToolUse hook on the Bash tool, through
cli.py's "guard" subcommand: the tool call arrives as JSON on stdin, exit
code 2 stops it, and stderr goes back into the session as the reason.
The rules are the board rules of the repository that vendors this
(its .claude/rules/board.md): the default branch takes only pull
requests, the board branch takes only the board workflow's writes, and
the one forced push the protocol allows is --force-with-lease on one's
own claude/ branch after a rebase; and the protocol's Claim rule, that the
claim branch of issue N is claude/N exactly, so a branch of other digits
after claude/, which reads as a claim and claims nothing, is not pushed.
A refusal is the rule working; the answer to one is a pull request, a
rebase or another branch, never a way around the hook.

The lease is also the one push GitHub itself would let take another
worker's claim: a worker that has fetched holds a remote-tracking ref
that matches the remote, so its --force-with-lease onto claude/N is
taken whoever claimed N, and a ruleset against forced pushes on claude/**
would refuse the prescribed lease too (spec/PROTOCOL.md, Work in the
open). So before a lease onto claude/N the guard reads whose the claim
is, from the claim file docs/agents/claims/N.md on the branch the lease
holds against (the commit a --force-with-lease=claude/N:<commit> names,
else, or when git cannot show that commit, as a lease written with a
shell variable names one, the remote-tracking ref refs/remotes/<remote>/
claude/N as the fetch left it), and whose claim the branch pushed makes
it, from the same file on that branch; on either side the file is read
at the tip, or, once the branch has deleted it as a finished branch has,
where its history last added it, the last so that an issue claimed again
after a merge is judged by its new claim and not by the one main's
history keeps. A lease whose two workers differ, or that would put a
branch with no claim file where another worker's claim stands, is
refused (other_holder). The file is read with git in the clone the push
runs in, so the fence costs a few git calls on a lease push and nothing
on any other command.

Each push is judged in each clone it may run in (refused_push,
push_places), since a session's shell is not always where its push goes:
git -C and env -C name another directory, and cd, pushd and popd move the
shell before a later command, so "cd ../other && git push origin HEAD"
pushes the other clone's branch, and its current and default branch, and
its claim files, are that clone's. The -C of git and env are joined on as
chdir reads them (push_place). The shell's moves are followed within one
command (placed_commands), each part of a pipeline and each list run with
& in a shell of its own, a cd or pushd to what is not a directory leaving
the shell where it stood, as the shell does when it fails there. But this
module cannot tell which of the moves it reads the shell makes: one in a
here-document's body, after a test that fails or in a function never
called counts as made, and one into a directory the command makes first,
as after mkdir, as failed. So a push after a cd, pushd or popd is judged
also where the command starts, as if no move ran, and, where a -C
relative to where the shell stands names its clone, in the session's own
clone too. Where this module cannot tell the directory, after cd - or
cd -P, with a path that holds what the shell expands, quoted or not,
since this module keeps no quotes, or that starts with bash's ~+, ~- or
~N (joined), or with the repository named by --git-dir, --bare, GIT_DIR
or GIT_REFERENCE_BACKEND, on the push or set before it, the push is
judged by what its words name and in the session's own clone. A push is
refused when any of these judgements refuses it. So a push is judged in
the session's own clone, as every push was judged before F21's fix, save
one whose clone git's or env's -C names by an absolute path, or a
relative one with no cd, pushd or popd before it in the same shell,
holding nothing the shell expands, and whose repository nothing names:
that push is judged in the clone the -C names alone. That is the clone it
runs in, save where the shell does what this module does not read, as
where the command first points a symbolic link on that path elsewhere
(ln -sfn), makes a repository there (git init, git worktree add) or
writes a .git file there, since the hook reads the path before the
command runs. This module does not trust a ~ to name a clone: it reads
one from the hook's own HOME, while a HOME the command sets before it
sends the shell elsewhere, and quotes, which this module does not keep,
leave a ~ as written, so a push by a -C from a ~ is judged where this
module reads it and in the session's own clone too (from_a_tilde). And it
does not follow a move or an export in a command handed to eval or trap
or in a file run by source or ., nor what a function or alias the
session's shell has before the command does, as one named git or cd may,
nor a variable that names the repository given its value by read,
printf -v or a for loop: it judges the push as if none of them ran. A
refusal reads as it did before the fix, and ends with what the hook could
not tell and how to have the push judged in its clone alone (UNMOVED,
UNTOLD) only where the reading that refused it is one the command's words
do not name.

The hook is the first fence and branch protection on GitHub the second,
so anything this module cannot read it lets through rather than wedging
the session: a push in quoting it misreads, a broken payload or a git that
answers nothing all pass, and GitHub still refuses what matters; a claim
file git cannot show, on either side of a lease, as on a branch whose
history never had one, is a claim the fence cannot judge, and the lease
passes as it did before the fence. Git is found behind what the shell
reads before a command: VAR=value assignments, the reserved words that
open a compound command or negate a pipeline (!, {, if, then, else, elif,
while, until, do, and function with its name), and the wrappers env,
command, nice, time, exec, nohup and timeout, each with its own options as
it reads them, a long one by any start of its name that names it alone, as
env --uns for env --unset, since GNU's getopt_long reads them so (WRAPPERS,
F26), env's -S string split into the words env splits it into, read
too with every ${NAME} in it unset, and after env's options, which a lone
- ends, the -, and every word that holds a =, which env reads as an
assignment where the shell would run it (F28 and its review), and git's
own options before the verb, each with its value where it takes one
(GIT_VALUE_OPTIONS), and a query after which git runs nothing
(GIT_QUERIES). A push is read as git reads its arguments: short options
in a bundle, as in -fu; a long one by any start of its name that names it
alone, as --force-w, or taken back with --no-; a value that is no
refspec; the refspecs that write many branches, the matching refspec :
and a pattern with a *; and each branch by the name git resolves,
refs/heads/main and heads/main both being main, and @ alone being HEAD.
git reads heads/X after the colon as X only where the remote has a branch
X, and makes a branch named heads/X where it has none. This module cannot
see the remote's branches and reads heads/X as X wherever it stands: that
is exact for the default branch and board, which are always there, and a
name the remote lacks becomes heads/X, a branch no rule names, so the
module at worst refuses a push git would have let make one, as
HEAD:heads/claude/03. What a command hands to a shell to run as text it
spells is read as a command of its own, from the shell that runs it
(handed_texts; F32): a shell's -c string, as in bash -c 'git push ...',
eval's words, trap's action, and a <<< here-string or what echo, printf
or cat prints into a pipe that a shell, or source /dev/stdin, reads. A
text it does not spell, from a file, a variable or a substitution, is a
value this module does not know and passes, as do a wrapper not named
here (sudo, xargs), a git alias for push, and a push git's configuration
sends where its words do not say. With push.default=upstream, or its old
name tracking, git 2.43 sends a branch named alone, as in "git push origin claude/7" and "git
push -u origin claude/7", and the current branch on a bare "git push", to
the branch it tracks: a claim branch made with "git switch -c claude/7
origin/main", which tracks main, went to main on each of the three, where
with --no-track the first push made claude/7. A remote's push refspecs
(remote.<name>.push) send a push elsewhere too. This module reads no
configuration.

shlex reads a newline as whitespace, which would glue a line onto the one
before it, so a small pass runs first (plain_text): it turns each unquoted
newline into the separator the shell reads it as, ends a comment at its
newline, and removes a backslash-newline, the shell's line continuation.
It tracks single quotes, double quotes and the backslash, and nothing else
of the shell: the quoting inside a command substitution is read as the
quoting around it, and bash's $'...' as a $ and a quoted span; so since
F24's fix each substitution's text is read apart as well (substituted,
below). bash runs
a command text one complete command at a time, so a quote that never
closes stops the line it opens on and what follows, and the commands
before it run; and it reads a here-document's body as text, whose quotes
it does not read. Until F25's fix this module read a command whole or not
at all: a quote that never closed hid every push in the command, those
bash ran before it among them, and an apostrophe in a body opened a quote
bash never opens, which one in a later body could close around a push, as
in the shape a session writes most, a commit message by here-document, a
push, then a pull request's body by another. So the hook reads a command
more than one way and judges every push any reading finds (readings,
placed_commands): whole, as before, with a quote that never closes read as
a character of its word, so the text on either side of it is read as
commands (plain_text, lenient), which reads too much of a line bash
refuses whole and never too little; and, where the command holds
here-documents this module tells (here_document_bodies), with their bodies
cut out, but for the command substitutions of a body the shell expands,
which it runs. The whole reading still reads each body as commands: a body
line that starts a push to main is refused, and the person rephrases it or
puts the text in a file, and a cd in a body moves the shell as this module
reads it, so a push after the body is judged also where the command
starts, as if no move ran (push_places), where the second reading judges
it too. Cutting bodies out in place of reading them was tried once and
dropped: telling a here-document from an arithmetic shift, a quoted marker
or an opener inside a substitution takes more of the shell than this hook
can read, and each misreading cut real commands out of sight. As a second
reading beside the whole, a body misread hides no command the whole
reading shows, and refuses too much at worst. The bodies told are those of
<< and <<- outside quotes and in a substitution, in double quotes or not,
as in "$(cat <<'EOF' ... )", each starting at the first newline at its
operator's level, or at the level around it where the substitution or
the double quotes it stands in close first, as in x=$(cat <<'EOF') (leave),
and ending past the line that is its word, a line of a body the shell
expands joined to the next where it ends in a backslash (continues),
where bash 5.2 ends one in a substitution (closes_substitution), or at
the end of the text; a <<< here-string, arithmetic, a comment and a
quoted span open none, a (( is arithmetic only where bash reads it so,
its parentheses counted past quotes (arithmetic_end), ((cd x); ...) being
a subshell in a subshell, and a case pattern's ) closes no substitution
(case_step), nor does one inside a parameter expansion, ${y:-)}, or a
backquoted span, which bash reads whole and whose here-documents lie
inside it (expansion_step, backquoted). What the shell runs of a body it
expands, its substitutions, is read as bash reads it, the bodies of its
own here-documents cut out in turn (substitutions_in), down to
MAX_NESTING, past which the command is refused unread (TOO_DEEP). Claude
Code reads a hook that runs out of time as no objection and runs the
command, so the reading is built to take a time in proportion to the
text: the scan matches each ( once (matching_parenthesis, scan), a cd
joins no path longer than the system takes (PATH_MAX) and a pushd stacks
in the same time at any depth (stacked); shlex takes a time in the square
of a word's length, so a command longer than MAX_COMMAND is refused
unread (TOO_LONG); and a command the guard still takes longer than
DEADLINE_SECONDS over, reading it and judging its pushes, git's answers
included, for a shape no test has timed, is refused too (TOO_SLOW),
never passed. The template's hook as it stands before it calls
refused_push (its tools/guard_push.py, which predates F28) calls
simple_commands, push_arguments and refusal itself, under neither
MAX_COMMAND nor DEADLINE_SECONDS: there a command too long or too slow to
read runs once the hook runs out of time, a gap older than F28 that goes
when the template calls refused_push. What F28 adds is bounded inside
push_arguments, so that hook meets it too: the -S strings by MAX_NESTING
and SPLIT_BUDGET, and a command it cannot read, which holds a push, comes
back as a word refusal refuses (Unreadable), never as an exception, which
that hook, exiting 1, would read as no objection (F28's second review).
A third reading, of each line alone, was weighed
and left out: with a quote that never closes read leniently it finds no
push in F25's shapes that the two miss, and what it alone would read is
text inside quotes that close, as a commit message given by -m "...",
which it would refuse wherever a line starts with a push; what hides a
push inside quotes that close is the quoting of a command substitution or
$'...', which want readings of their own. A backquoted span is read
whole in the lenient readings, as one part of its word (backquoted_word),
so a quote inside one pairs with nothing outside it (F25's fourth
review). And bash runs the text of every command substitution, $( ) or a
backquoted span, wherever it expands one, in double quotes, a parameter
expansion and arithmetic too, with quoting of its own, and the text of a
process substitution, <( ) or >( ), outside double quotes, where shlex
reads <( as a redirection that takes the next word: until F24's fix
out="$(git push ... 2>&1)" and cat <(git push ...)
passed. So the text of each substitution bash runs is read apart, as the
commands of a subshell after the reading, those nested in it inside the
same subshell after its own text, a backquoted span's with the
backslashes bash takes off taken off (substitutions, substituted), from
the command with its here-documents' bodies cut out, so a quoted body in
a substitution, a commit message's, stays text, as F25's fix made it. A
substitution is judged where the shell stands after the whole reading,
and, after a move, where the command starts, as any push after a move
is, so a cd after it in the command may judge it in a clone it does not
run in, which refuses too much at worst. The scans, of here-documents
and of substitutions, step over a $'...' whole with its escapes where
bash reads commands and in a ${ } outside double quotes (ansi_c_end), so
one hides nothing after it; read a single quote in arithmetic, and in a
${ } in double quotes after an operator other than a pattern's, as bash
does, pairing with the next but quoting nothing (APOSTROPHE, brace_state);
and read arithmetic to its closing parentheses in a context of its own
(ARITHMETIC), each since F24's first review.

And bash builds words: it decodes $'...' and reads $"..." as double
quotes, runs a substitution inside a word and splices what it prints
there, drops an unquoted word that expands to nothing, so $NOPE git push
runs git, and expands braces and globs; until F29's fix the guard read
each word as written, so git $'pu\x73h', git pu$(true)sh, $NOPE git
push, git {push,} and git push origin ma$(echo i)n passed. Since then
each $'...' is read as the text it decodes to (ansi_c_text), the command
is read once more as its own text with each substitution a hole in its
word (read_apart), and each simple command's words are read as written,
with their braces expanded, and with every variable unset and every
substitution printing nothing (word_ways); a word with holes is read as
any text it may be built into (may_be), so a first word that may be git
and a verb that may be push are read as git and push, an option whose
name bash builds is refused, since it may be --force, and so is a branch
name bash may build into the default or the board branch; a word that
is one expansion whole is still read as a value the guard cannot know.
Since its first review so is one made of expansions alone, as
"${PREFIX}${N}", and a refspec of them, as $B:$B; a built option passes
after a whole flag no other extends, as --no-verify${NOPE}; the hook's
gate takes line continuations out; coproc runs the command after it; an
unquoted ${ } a blank stands in goes on into the words after it; a ${ }
reads an array's subscript and decodes a $'...' in its default; and an
extended pattern, as @(main), is read as part of its word too.
A GIT_DIR=x alone, which git sees only where it was exported before, may
name the repository, so the push after it is judged both ways. The scan
of here-documents keeps a word's start across a line continuation, so a
# after one opens a comment, steps into arithmetic, whose substitutions
bash runs, and reads a here-document's delimiter word whole, a $( ) in it
text, as bash does.

What this leaves: the whole reading still reads each body as commands
with an ordinary line's quoting, so a body line that holds a double quote
before a substitution naming a push may be refused, as before F24's fix,
and so may a body line that names a push after an expansion, as one that
starts with the push is;
a substitution's own text is made plain as any text is, so a push named
inside a $'...' that holds an escaped apostrophe, in a substitution, may
be refused; a push named in a here-document's delimiter word, which bash
reads as text, is still read by the whole reading and refused, and so is
one in a ${ } whose word bash never expands or cannot parse; a quote the
guard drops, as in "${NOPE}" git push, which bash runs as a command named
nothing, reads as an expansion that may give nothing; each refuses too
much. A substitution or variable that is a word whole, as git push origin
$(echo main), is read as a value the guard cannot know and passes, as
$BRANCH does; a hole is never read as giving more than one word, a : or a
+. Branch protection stands behind each.
"""

import bisect
import contextvars
import fnmatch
import functools
import itertools
import os
import re
import shlex
import subprocess
import time
import unicodedata
import weakref
from typing import NamedTuple

# git push's own options, as its parse-options reads them (git 2.43's
# git push -h): each long one by what it takes, VALUE a value attached with
# = or given as the next word, which is then no refspec, OPTIONAL a value
# only when attached with =, FLAG none. All but NOT_NEGATED also read as
# --no-<name>, which takes the option back, and git takes any start of a
# long name that names one option alone for the whole of it, --force-w for
# --force-with-lease and --mi for --mirror. The short ones bundle, -fu
# being -f and -u, and -o takes the rest of its word, or the next word.
FLAG, VALUE, OPTIONAL = "flag", "value", "optional"
PUSH_OPTIONS = {
    "verbose": FLAG, "quiet": FLAG, "repo": VALUE, "all": FLAG, "branches": FLAG, "mirror": FLAG,
    "delete": FLAG, "tags": FLAG, "dry-run": FLAG, "porcelain": FLAG, "force": FLAG,
    "force-with-lease": OPTIONAL, "force-if-includes": FLAG, "recurse-submodules": VALUE, "thin": FLAG,
    "receive-pack": VALUE, "exec": VALUE, "set-upstream": FLAG, "progress": FLAG, "prune": FLAG,
    "no-verify": FLAG, "verify": FLAG, "follow-tags": FLAG, "signed": OPTIONAL, "atomic": FLAG,
    "push-option": VALUE, "ipv4": FLAG, "ipv6": FLAG,
}
NOT_NEGATED = {"no-verify", "verify", "ipv4", "ipv6"}
# Every whole long name git push reads, by whether it may take a value:
# each option's own, and no-<option> for each one git lets be negated,
# which takes none (whole_flag).
FLAG_NAMES = frozenset({name for name, kind in PUSH_OPTIONS.items() if kind == FLAG}
                       | {f"no-{name}" for name in PUSH_OPTIONS if name not in NOT_NEGATED})
VALUE_NAMES = frozenset(name for name, kind in PUSH_OPTIONS.items() if kind != FLAG)
PUSH_LETTERS = {"v": "verbose", "q": "quiet", "d": "delete", "n": "dry-run", "f": "force",
                "u": "set-upstream", "o": "push-option", "4": "ipv4", "6": "ipv6"}
# The options that push branches wholesale rather than the one named.
WHOLESALE = ("all", "branches", "mirror")

# Global options of git itself that take a separate value, between "git"
# and the verb, as in "git -C clone push": git 2.43's (git.c,
# handle_options), --config-env, --attr-source and --shallow-file among
# them, whose value was read as the verb until F21's fix and its review,
# so "git --config-env user.name=HOME push origin main" read as no push.
GIT_VALUE_OPTIONS = {"-C", "-c", "--config-env", "--git-dir", "--work-tree", "--namespace", "--attr-source",
                     "--shallow-file"}
# Global options after which git prints what it is asked and exits, running
# no command: "git --exec-path push origin main" pushes nothing. So does
# --list-cmds=<group>.
GIT_QUERIES = {"--exec-path", "--html-path", "--man-path", "--info-path"}
# git 2.43's other global options (git.c, handle_options), with no value or
# with one attached by =: any other option before the verb git refuses, and
# runs nothing. verb_at reads them only after a first word that bash may
# build into git and that is not git as written (F29), since a value the
# guard cannot know stands as git there only in a world where it is git.
GIT_FLAGS = {"-p", "--paginate", "-P", "--no-pager", "--no-replace-objects", "--bare", "--literal-pathspecs",
             "--no-literal-pathspecs", "--glob-pathspecs", "--noglob-pathspecs", "--icase-pathspecs",
             "--no-optional-locks", "--no-lazy-fetch", "--no-advice"}
GIT_ATTACHED = ("--exec-path=", "--git-dir=", "--namespace=", "--work-tree=", "--super-prefix=", "--config-env=",
                "--list-cmds=", "--attr-source=")
# How a push names its repository by a path, an option of git's own,
# attached or not, or a variable, set on the push, before it alone or
# exported earlier: the push then goes where this module does not look
# (push_place, Shell.named_repository). --work-tree names a work tree and
# leaves the repository to the directory. GIT_REFERENCE_BACKEND, which git
# reads since 2.54, keeps the refs where its files://<path> says, so a
# lease holds against that path's refs wherever -C points (F21's fourth
# review: git 2.55 took a lease by git -C into a clone that had not fetched
# the claim, against the refs of one that had); GIT_COMMON_DIR moves no
# ref, since git finds its ref store without reading it (git 2.43 and 2.55
# both refused that lease as stale), and is not counted.
GIT_REPOSITORY_OPTIONS = {"--git-dir", "--bare"}
GIT_REPOSITORY_VARIABLES = ("GIT_DIR", "GIT_REFERENCE_BACKEND")
# What a substitution stands as in the text of the one around it, once its
# own text is read apart (substitutions): a character of Unicode's private
# use area, which no command a session writes holds, so the reading of a
# word tells where a substitution stood in it and reads it as a hole bash
# fills (word_pieces; F29: it stood as $x until then, so pu$(true)sh read
# as pu and a variable named xsh). A path or a branch it gives is one this
# module cannot tell, as it was before F24's fix.
SUBSTITUTED = "\ue000"
# What the shell expands in a path before cd or git sees it: a variable, a
# command's output, a glob; such a path names a directory this module does
# not know (joined).
EXPANDED = set("$`*?[" + SUBSTITUTED)
# bash's tilde prefixes that name the shell's own directories, followed by
# a / or nothing: ~+ the directory it stands in, ~- the one it stood in
# before, and ~N, ~+N and ~-N an entry of its directory stack. A path that
# starts with one names a directory this module does not know (joined),
# where os.path.expanduser left it a directory named ~+ or ~- (F21's
# fourth review).
SHELL_TILDE = re.compile(r"~(?:[+-][0-9]*|[0-9]+)(?:/|$)")
# The longest path the system takes in one call, in bytes with its NUL
# (Linux's PATH_MAX): a longer one it refuses (ENAMETOOLONG), where the
# shell's cd still reaches it a step at a time, so this module cannot tell
# the directory such a path names (joined). The bound also keeps a long
# run of cd from joining ever longer paths, a time in the square of the
# run (F25's third review: 43,000 of "cd a;" took 16 seconds).
PATH_MAX = 4096


class Wrapper(NamedTuple):
    """How a command that runs another reads its own words first, as GNU's
    getopt_long does: its short options that take a value, attached as in
    -n5 or as the next word; every long option it has, each by what it
    takes, as git push's are read (FLAG, VALUE, OPTIONAL), the flags among
    them, since a word names an option by any start of its name that names
    no other (long_option_of); and how many words follow its options before
    the command it runs, as timeout's duration does."""

    short_values: str
    long_options: dict
    operands: int = 0


# The two long options every GNU program reads.
GNU_STANDARD = {"help": FLAG, "version": FLAG}
# The commands a person puts before git to run it another way, with every
# option each reads before the command (F26: the guard compared a long one
# with the names that take a value, in full, so "env --uns HOME", "nice
# --adj 5" and "timeout --sig KILL" were read as running HOME, 5 and the
# duration). GNU coreutils' env, nice, nohup and timeout as 9.4's --help
# lists them, with env's -a and --argv0 of 9.5, its --env0-from of 9.12
# and its --quoting-style, which coreutils' NEWS lists for the release
# after 9.12, since a version that lacks an option refuses the word and
# runs nothing, and one that has it runs the command after it; GNU time as
# its src/time.c reads them in 1.9 and 1.10, whose manual's --output is a
# start of --output-file; and bash's command and exec, which take no long
# option (bash's own time takes only -p, a flag). env's short options that
# take a value hold the BSD ones too: -P altpath, the directories the env
# of FreeBSD and macOS searches for the command, and FreeBSD's -L and -U,
# a login class and a user, since those makes run the command after the
# value; GNU's env refuses all three and runs nothing (exit 125), so
# reading the next word as their value costs only the refusal of a push
# GNU would not have run. Any other option another make's wrapper reads
# and GNU's does not is read as a flag, and the push behind it passes, as
# before F26's fix.
WRAPPERS = {
    "env": Wrapper("uCSaPLU", {"argv0": VALUE, "ignore-environment": FLAG, "env0-from": VALUE, "null": FLAG,
                            "unset": VALUE, "chdir": VALUE, "default-signal": OPTIONAL, "ignore-signal": OPTIONAL,
                            "block-signal": OPTIONAL, "list-signal-handling": FLAG, "debug": FLAG,
                            "quoting-style": VALUE, "split-string": VALUE, **GNU_STANDARD}),
    "command": Wrapper("", {}),
    "nice": Wrapper("n", {"adjustment": VALUE, **GNU_STANDARD}),
    "time": Wrapper("of", {"append": FLAG, "format": VALUE, "output-file": VALUE, "portability": FLAG, "quiet": FLAG,
                           "verbose": FLAG, **GNU_STANDARD}),
    "exec": Wrapper("a", {}),
    "nohup": Wrapper("", GNU_STANDARD),
    "timeout": Wrapper("ks", {"foreground": FLAG, "kill-after": VALUE, "preserve-status": FLAG, "signal": VALUE,
                              "verbose": FLAG, **GNU_STANDARD}, operands=1),
}
# command -v and -V say what would run and run nothing.
DESCRIBES_ONLY = set("vV")
# A word the shell reads as an assignment to the command's environment,
# VAR=value or VAR+=value.
ASSIGNMENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\+?=")
# The reserved words the shell reads before a command without running
# anything of their own: those that open or go on with a compound command,
# and ! that negates a pipeline. The command after one is a command of its
# own, so "if git push ...", "while ! git push ...", "do git push ..." and
# "{ git push ...; }" each run the push. function takes the next word as
# the name it defines, and coproc, which runs the command after it, takes
# it as the coprocess's name before a compound command (coproc_named).
RESERVED = {"!", "{", "if", "then", "else", "elif", "while", "until", "do"}

# A claim branch is claude/ and the issue's number as GitHub writes it,
# one to nine ASCII digits with no leading zero, the rule claims.NUMBER
# states (this module imports nothing of the package, since the plugin
# marketplace takes it alone). Any other digits after claude/, as in
# claude/03, read as a claim and claim nothing, since GitHub's lock, the
# refusal of a second push, works per name.
CLAIM_NUMBER = re.compile(r"0|[1-9][0-9]{0,8}")
# The claim file the protocol's Claim rule names, docs/agents/claims/N.md,
# and the line in it that names the claim's worker, "- Worker: name" (or
# the older "- Account: name"), with a note in parentheses after the name
# left to the reader, as claims.WORKER_LINE and claims.TRAILING_NOTE read
# it; tests/test_guard.py holds the two to each other.
CLAIM_DIR = "docs/agents/claims"
WORKER_LINE = re.compile(r"^\W*(?:worker|account)\W*:\s*(.+?)\s*$", re.IGNORECASE)
TRAILING_NOTE = re.compile(r"\s*\(.*\)\s*$")
# Where a fetch leaves a remote's branches in the clone, which a lease
# that names no commit holds against; origin unless the push names a
# remote (remote_ref).
REMOTE_REFS = "refs/remotes/"

# What ends a word for the shell: a # right after one of these starts a
# comment; anywhere else it is part of the word. A here-document's
# delimiter word ends at one too.
WORD_BREAKS = " \t\n;&|()<>"
# A double-quoted span as the shell reads it, up to the first " that no
# backslash escapes, a backslash and the character after it being one pair
# wherever they stand; and such a pair.
DOUBLE_QUOTED = re.compile(r'"(?:[^"\\]|\\.)*"', re.DOTALL)
BACKSLASH_PAIR = re.compile(r"\\(.)", re.DOTALL)
# The pairs whose backslash quote removal takes off inside double quotes,
# as bash reads a here-document's delimiter written there; a backslash
# before anything else is kept.
DOUBLE_QUOTED_ESCAPE = re.compile(r'\\([$`"\\\n])')

# The longest command the hook reads, in characters: shlex builds each
# word a character at a time, in a time that grows with the square of the
# word's length (a word of 800,000 characters took 9 seconds), so a longer
# command is refused unread rather than run the hook out of time, which
# would pass it.
MAX_COMMAND = 262144
TOO_LONG = (f"a command of more than {MAX_COMMAND} characters is refused unread: the guard would take longer to "
            "read it than a hook may run; write the text to a file first (.claude/rules/board.md)")
# How deep the here-documents of a command may nest, each in a command
# substitution of a body the shell expands, before the hook refuses the
# command unread: bash reads each such body again inside the one around
# it, and so would this module, the time it takes growing with the depth.
MAX_NESTING = 16
# The refusal of a command whose here-documents nest deeper than that,
# which the hook refuses unread rather than run out of time and pass it.
TOO_DEEP = (f"a command whose here-documents nest in command substitutions more than {MAX_NESTING} deep is "
            "refused: the guard cannot read it to tell which of its lines run; write the text to a file "
            "first (.claude/rules/board.md)")

# How long the guard may take over one command, reading it and judging
# its pushes, git's answers included, in seconds of the monotonic clock,
# before the hook refuses it: Claude Code reads a hook that runs out of
# time as no objection and runs the command, so a guard that runs long
# must end in a refusal, never in the hook's timeout (60 seconds unless a
# settings file says otherwise). Each part of the reading is meant to take
# a time in proportion to the text, and a command at MAX_COMMAND reads in
# about a second; the deadline stands behind that for a shape no test has
# timed (F25's third review found two that took minutes), and bounds the
# git calls the judging makes, a few for each clone a push may run in, of
# which a command can name thousands (its fourth review: 11,000 pushes,
# each by git -C into a clone of its own, took 42 seconds).
DEADLINE_SECONDS = 10
TOO_SLOW = (f"a command the guard cannot read and judge in {DEADLINE_SECONDS} seconds, or for which git did not "
            "answer in time, is refused: a hook that runs out of time lets the command run; write the text to a file "
            "first, or push in a command of its own (.claude/rules/board.md)")

UNREADABLE = ("a command the guard cannot read, for the size of what it builds, is refused: a hook that fails lets "
              "the command run; write the text to a file first, or push in a command of its own "
              "(.claude/rules/board.md)")

# The hook's exit codes: 0 lets the tool call run, 2 stops it and sends
# stderr back into the session.
PASS = 0
REFUSE = 2


def quote_end(command, i):
    """Where the quoted span that opens at command[i], a ' or a ", ends,
    past its closing quote: a single quote's at the next ', a double
    quote's at the first " that no backslash escapes; None when it never
    closes."""
    if command[i] == "'":
        end = command.find("'", i + 1)
        return None if end == -1 else end + 1
    match = DOUBLE_QUOTED.match(command, i)
    return match.end() if match else None


def as_quoted(span):
    """A quoted span as shlex is handed it: a double-quoted one without its
    backslash-newlines, the line continuation the shell removes there too,
    and with \\$ and \\` made $ and `, as bash makes them there and shlex
    does not (F28's second review: "env -S \"\\${NOPE}git push origin
    main\"" hands env ${NOPE}, which the guard read as \\${NOPE}), every
    other pair kept; a single-quoted one as written."""
    if not span.startswith('"'):
        return span
    return BACKSLASH_PAIR.sub(lambda pair: QUOTED_PAIRS.get(pair.group(1), pair.group(0)), span)


# What a backslash pair inside double quotes is handed to shlex as, where
# it is not kept whole (as_quoted).
QUOTED_PAIRS = {"\n": "", "$": "$", "`": "`"}

# bash's ANSI-C quoting, $'...', as its ansicstr decodes the text (F29):
# the escapes of one letter, \nnn in octal, \xHH or \x{...} in hex, \uHHHH
# and \UHHHHHHHH as a character, and \cX as a control character; any other
# backslash is kept with the character after it, and the text ends at a NUL
# byte, as a C string does.
ANSI_C_ESCAPES = {"a": "\a", "b": "\b", "e": "\x1b", "E": "\x1b", "f": "\f", "n": "\n", "r": "\r", "t": "\t",
                  "v": "\v", "\\": "\\", "'": "'", '"': '"', "?": "?"}
OCTAL = re.compile(r"[0-7]{1,3}")
HEX_DIGITS = {"x": re.compile(r"[0-9A-Fa-f]{1,2}"), "u": re.compile(r"[0-9A-Fa-f]{1,4}"),
              "U": re.compile(r"[0-9A-Fa-f]{1,8}")}
BRACED_HEX = re.compile(r"\{([0-9A-Fa-f]*)\}")


def ansi_c_close(command, i):
    """Where the $'...' whose text starts at command[i], past its $', ends:
    past the apostrophe no backslash escapes, since in it a backslash
    escapes the character after it, an apostrophe too; None where none
    closes it."""
    while i < len(command):
        if command[i] == "\\":
            i += 2
        elif command[i] == "'":
            return i + 1
        else:
            i += 1
    return None


def ansi_c_text(text):
    """What bash's $'...' gives for the text between its quotes (ANSI_C_ESCAPES):
    git $'pu\\x73h' origin $'ma\\151n' runs git push origin main (F29).
    A NUL ends it, as it ends any word bash builds."""
    return escapes_read(text).partition("\0")[0]


def escapes_read(text):
    """The text with its backslash escapes read as $'...' reads them, a NUL
    kept where one is written: echo and printf print it, and a shell that
    reads their output from a pipe drops it and reads on (printed_texts;
    F32's second review)."""
    out, i = [], 0
    while i < len(text):
        ch = text[i]
        if ch != "\\" or i + 1 == len(text):
            out.append(ch)
            i += 1
            continue
        escape, i = text[i + 1], i + 2
        if escape in ANSI_C_ESCAPES:
            out.append(ANSI_C_ESCAPES[escape])
        elif escape in "01234567":
            digits = OCTAL.match(text, i - 1)
            out.append(chr(int(digits.group(0), 8) & 0xFF))
            i = digits.end()
        elif escape == "x" and (braced := BRACED_HEX.match(text, i)):
            out.append(chr(int(braced.group(1) or "0", 16) & 0xFF))
            i = braced.end()
        elif escape in HEX_DIGITS and (digits := HEX_DIGITS[escape].match(text, i)):
            out.append(code_point(int(digits.group(0), 16), escape, digits.group(0)))
            i = digits.end()
        elif escape == "c" and i < len(text):
            control, upper = text[i], text[i].upper()
            i += 2 if text.startswith("\\\\", i) else 1
            out.append("\x7f" if control == "?" else chr(ord(upper if len(upper) == 1 else control) & 0x1F))
        else:
            out.append("\\" + escape)
    return "".join(out)


def code_point(value, escape, digits):
    """The character a \\x, \\u or \\U escape gives, its escape as written
    where the value names no character."""
    try:
        return chr(value)
    except (ValueError, OverflowError):
        return "\\" + escape + digits


def ansi_c_word(command, i):
    """The $'...' that opens at command[i] as one single-quoted span of the
    text it decodes to (ansi_c_text), and where it ends; None where none
    opens there or none closes."""
    if not command.startswith("$'", i):
        return None
    end = ansi_c_close(command, i + 2)
    return (single_quoted(ansi_c_text(command[i + 2:end - 1])), end) if end else None


def locale_quoted(command, i):
    """Whether a $"..." opens at command[i]: bash reads it as the double-quoted
    span after the $, where no message catalog translates it."""
    return command.startswith('$"', i)


def joins_lines(pair):
    """Whether a backslash pair is a line continuation, which bash takes out
    before it reads on, so whether a word starts is as it was before it."""
    return pair == "\\\n"


def single_quoted(text):
    """text as one single-quoted word for shlex, an apostrophe in it closed,
    escaped and opened again."""
    return "'" + text.replace("'", "'\"'\"'") + "'"


def plain_text(command, lenient=False, spans_run=True, patterns=False):
    """The command with each unquoted newline made a ; so that shlex, which
    reads a newline as whitespace, ends a command there as the shell does.
    A comment ends at its newline, which then separates; a backslash-newline
    is the line continuation the shell removes; a newline inside quotes
    stays inside its word. A quote that never closes leaves the rest of the
    text as written, and a backslash that ends it stays alone, so shlex
    refuses the whole (plain_shell); read leniently, each is a character of
    its word, as if a backslash escaped it, and the text after the quote is
    read as text no quote opened, so the commands on either side of it are
    read (F25: bash runs each complete command before the line it cannot
    read, and that line's own as it reads them). A lenient reading reads a
    backquoted span as one part of its word (backquoted_word) and, where
    spans_run, its text after the rest as the commands of a subshell.

    Since F29's fix bash's quotings of its own are read as it reads them: a
    $'...' as the text it decodes to (ansi_c_text), one single-quoted span,
    where until then it read as a $ and a quoted span, so git push origin
    $'main' read as a push to $main; a $"..." as the double-quoted span it
    is, where no message catalog translates it; $$, the shell's process id,
    as a parameter whole, so the quote after it opens no $'...'; and a line
    continuation leaves whether a word starts where it was, since bash takes
    it out before it reads on, so a # after one that follows a blank opens a
    comment (F29: echo a \\, a newline and #' hid the push on the next line
    inside a quote bash never opened). Where patterns, each extended
    pattern, @( ), ?( ), *( ), +( ) or !( ), is kept whole as a part of its
    word (pattern_group), as bash reads one with extglob on, where without
    it the parentheses end the word (pattern_reading)."""
    out, spans, word_start, i = [], [], True, 0
    closes = closing_parentheses(command) if patterns else None
    while i < len(command):
        ch = command[i]
        if ch == "$" and command.startswith("$$", i):
            out.append("$$")
            i, word_start = i + 2, False
        elif ch == "$" and (ansi_c := ansi_c_word(command, i)):
            out.append(ansi_c[0])
            i, word_start = ansi_c[1], False
        elif ch == "$" and locale_quoted(command, i):
            i, word_start = i + 1, False
        elif ch in "'\"":
            end = quote_end(command, i)
            if end is None and not lenient:
                out.append(as_quoted(command[i:]))
                break
            out.append("\\" + ch if end is None else as_quoted(command[i:end]))
            if lenient and spans_run and end and ch == '"':
                spans += [plain_text(inner, lenient=True) for inner in backquoted_in(command[i + 1:end - 1])
                          if inner.strip()]
            i, word_start = end or i + 1, False
        elif ch == "\\":
            pair = command[i:i + 2]
            out.append("" if pair == "\\\n" else "\\\\" if pair == "\\" and lenient else pair)
            i, word_start = i + 2, word_start and joins_lines(pair)
        elif ch == "#" and word_start:
            end = command.find("\n", i)
            i = len(command) if end == -1 else end
        elif patterns and (group := pattern_group(command, i, closes)):
            out.append(group[0])
            i, word_start = group[1], False
        elif ch == "`" and lenient and (span := backquoted_span(command, i)):
            out.append(backquoted_word(span.group(0)))
            if spans_run and span.group(0)[1:-1].strip():
                spans.append(plain_text(span.group(0)[1:-1], lenient=True))
            i, word_start = span.end(), False
        else:
            out.append(" ; " if ch == "\n" else ch)
            i, word_start = i + 1, ch in WORD_BREAKS
    return "".join(out) + "".join(f" ; ( {inner} ) " for inner in spans)


def pattern_group(command, i, closes):
    """The extended pattern that opens at command[i], one of @?*+! and a
    parenthesis not after a $, as the text shlex reads as one part of a
    word, each character that would end it escaped, and where it ends;
    None where none opens there, or its parenthesis does not close on its
    line. closes is the command's parentheses' (closing_parentheses)."""
    if command[i] not in "@?*+!" or not command.startswith("(", i + 1) or command[i - 1:i] == "$":
        return None
    end = closes.get(i + 1)
    if end is None:
        return None
    return command[i] + "".join("\\" + ch if ch in ENDS_A_WORD else ch for ch in command[i + 1:end]), end


def closing_parentheses(text, escapes=True):
    """Where the parenthesis that closes each ( of text ends, by the index
    of the (, read in one pass with a stack, so a run of ( that never close
    costs no more than one: none closes past a newline, and, where escapes,
    a backslash takes the character after it out."""
    closes, stack, i = {}, [], 0
    while i < len(text):
        ch = text[i]
        if ch == "\\" and escapes:
            i += 2
            continue
        if ch == "(":
            stack.append(i)
        elif ch == ")" and stack:
            closes[stack.pop()] = i + 1
        elif ch == "\n":
            stack.clear()
        i += 1
    return closes


# What ends a word for shlex outside quotes, which an extended pattern keeps
# (pattern_group).
ENDS_A_WORD = frozenset(" \t()|&;<>")


def backquoted_span(command, i):
    """The match of the backquoted span that opens at command[i]
    (BACKQUOTED), or None where none opens there or it never closes."""
    return BACKQUOTED.match(command, i) if command[i] == "`" else None


def backquoted_in(text):
    """The text of each backquoted span in text, a double-quoted span's
    inside, past backslash pairs and past each $( ) whole: bash runs each
    such span when it expands the double quotes, so a lenient reading reads
    each as commands too (plain_text), but what stands in a $( ) is that
    substitution's own text, where a backquote may sit in single quotes or
    in a quoted here-document's body, which bash never runs (F25's fifth
    review: "$(cat <<'EOF' ... `git push --force` ... EOF)", a commit
    message naming a command, was refused), so the substitution is read as
    0.5.0 reads it, one part of its word. Single quotes inside double
    quotes are text, so a span between them is found and runs."""
    found, memo, i = [], {}, 0
    while i < len(text):
        if text.startswith("$(", i):
            close = matching_parenthesis(text, i + 1, memo)
            if close is None:
                break
            i = close + 1
            continue
        span = BACKQUOTED.match(text, i) if text[i] == "`" else None
        if span:
            found.append(span.group(0)[1:-1])
        i = span.end() if span else i + (2 if text[i] == "\\" else 1)
    return found


def backquoted_word(span):
    """A backquoted span, ` to the first ` no backslash escapes
    (BACKQUOTED), as a lenient reading hands shlex it in place: one part of
    a word, single-quoted, with nothing inside it read, since bash reads
    the span whole and a quote inside it pairs with nothing outside (F25's
    fourth review: x=`echo "`; git push ... # " read the push inside a
    quoted word). The reading reads the span's text as well, as the
    commands of a subshell after the rest (plain_text), so a push inside
    it, which bash runs, is judged too, while the word it stands in, as
    git -C `pwd` push, stays whole."""
    return "'" + span.replace("'", "'\"'\"'") + "'"


def plain_shell(command):
    """The command made plain for shlex (plain_text), a quote that never
    closes left for shlex to refuse whole, as shell_parts reads it."""
    return plain_text(command)


# --- here-documents ------------------------------------------------------------

class HereDocument(NamedTuple):
    """A here-document's operator as the shell reads it: the word whose line
    ends its body, its quotes taken off; whether the body's lines lose their
    leading tabs, as after <<-; whether the shell expands the body, as it
    does where no part of the word is quoted; and the level the operator
    stands at (Context), where the first newline the shell reads after it
    starts its body, None at the top."""

    delimiter: str
    strip_tabs: bool
    expands: bool
    within: object


class Body(NamedTuple):
    """A here-document's body, as the span of the command from start to
    end, and what of it the shell runs: the command substitutions of a body
    it expands (substitutions_in)."""

    start: int
    end: int
    runs: tuple = ()


class TooDeep(Exception):
    """A command whose here-documents nest in the substitutions of bodies
    the shell expands deeper than this module reads (MAX_NESTING)."""


class TooSlow(Exception):
    """A command the guard ran past its deadline over (DEADLINE_SECONDS)."""


# The monotonic clock the deadline is read by, and the deadline of the
# guard's run under way, None where none is (refused_push sets it).
clock = time.monotonic
DEADLINE = contextvars.ContextVar("deadline", default=None)


def seconds_left():
    """The seconds left before the deadline of the run under way, None
    where none is under way."""
    deadline = DEADLINE.get()
    return None if deadline is None else deadline - clock()


def in_time():
    """Raise TooSlow where the run under way has passed its deadline;
    called at the head of each loop of the reading (scan, parts_of,
    placed_parts), before the judging asks about each clone, and before
    each git call (git_line), so no shape of command holds the hook past
    it."""
    left = seconds_left()
    if left is not None and left < 0:
        raise TooSlow()


# What the scan of a command's here-documents tells levels by: double
# quotes, a command substitution, a parameter expansion, ${ }, and
# arithmetic, inside which a newline starts no body, since bash reads the
# expansion whole as part of a word before it reads the newline after it.
LEVELS = ('"', "$(", "${", "$[", "$((")
# A backquoted span as bash's parse_matched_pair reads it, up to the first
# ` that no backslash escapes, with nothing nested inside (F25's third
# review: a case pattern's ) inside one closed the substitution around it).
BACKQUOTED = re.compile(r"`(?:[^`\\]|\\.)*`", re.DOTALL)
# The contexts inside which a << is text: double quotes, and a parameter
# expansion or old-style arithmetic, $[ ], with a [ nested in it, each of
# which bash reads whole as part of a word, and inside those an apostrophe
# that pairs but quotes nothing (APOSTROPHE), and arithmetic, $(( )) or
# (( )), read apart (substitutions).
APOSTROPHE = "'"
ARITHMETIC = "$(("
TEXT_ONLY = ('"', "${", "$[", "[", APOSTROPHE, ARITHMETIC)
# The characters of a parameter expansion's operator, after its name, as
# bash's parse_matched_pair tells its state (brace_state): after one that
# removes or replaces a pattern, #, %, /, ^ or , after the name, a single
# quote quotes inside double quotes too; after any other, :- among them,
# it pairs with the next but quotes nothing there, so bash 5.2 runs the
# push in echo "${v:-'$(git push ...)'}" (F24's first review).
BRACE_OPERATORS = "#%^,~:-=?+/"
PATTERN_OPERATORS = "#%^,/"
# The reserved words after which a word is still a command's first, so a
# case there is a case command.
BEFORE_A_COMMAND = ("then", "do", "else", "elif", "if", "while", "until", "!", "{", "time")
# The characters those words open with, so the scan looks for one only
# where a word may be one: a run of parentheses, each a word's start, was
# read against every word (F29's second review: @ and 262,000 of them took
# 3.5 seconds).
BEFORE_A_COMMAND_OPENS = frozenset(word[0] for word in BEFORE_A_COMMAND)


class Context:
    """One context the scan of a command's here-documents stands in, inside
    the one around it (outer): double quotes ("), a command substitution
    ($(), a parenthesis inside one ((), which the scan counts to tell the
    ) that closes the substitution, or a case command inside one, before
    its in (case), among its patterns (in) or in a pattern's commands
    (case)), since bash reads a pattern's ) as closing no substitution
    (F25), or a parameter expansion (${), inside which a ( or ) is text
    (F25's third review). Its level is the innermost double quotes,
    substitution or parameter expansion at or around it, where a
    here-document's body starts at a newline, and whether that stands in a
    substitution. A context is one object, told from another by identity,
    so opening or closing one takes the same time however deep the scan
    stands."""

    __slots__ = ("kind", "outer", "level", "in_substitution", "quoted", "state", "seen", "end")

    def __init__(self, kind, outer=None, at=0, end=None):
        self.kind, self.outer = kind, outer
        self.level = self if kind in LEVELS else level(outer)
        self.in_substitution = kind == "$(" or bool(outer and outer.in_substitution)
        self.quoted = kind in ("${", APOSTROPHE) and outer is not None and (
            outer.kind == '"' or (outer.kind in ("${", APOSTROPHE) and outer.quoted))
        self.state, self.seen, self.end = "first", at, end

    def __iter__(self):
        """The kinds of this context and of each around it, outward."""
        context = self
        while context is not None:
            yield context.kind
            context = context.outer


def level(contexts):
    """The double quotes or substitution the scan stands in at contexts
    (Context), None at the top: a parenthesis or a case inside a
    substitution opens no level of its own, as a subshell's does not at the
    top."""
    return None if contexts is None else contexts.level


def takes_here_document(contexts):
    """Whether a << that the scan of here_document_bodies meets in contexts,
    the double quotes, substitutions and parentheses around it, opens a
    here-document: outside quotes, and in a substitution inside double
    quotes, as in "$(cat <<'EOF' ...)", the form a session writes a commit
    message or a pull request's body in; in double quotes themselves it is
    text, and so it is in a parameter expansion, ${x:-<<EOF}, which bash
    reads whole as part of a word."""
    return contexts is None or contexts.kind not in TEXT_ONLY


def backquote_end(command, i):
    """Where the backquoted span that opens at command[i] ends, past the `
    that closes it (BACKQUOTED), or the end of the text where none does."""
    match = BACKQUOTED.match(command, i)
    return match.end() if match else len(command)


def backquoted(command, i, depth):
    """Where the scan goes on past the backquoted span at command[i], and
    the bodies of the here-documents inside it, or None where command[i]
    opens none. bash reads the span whole, to the first ` no backslash
    escapes, so a ) or a newline inside it closes and starts nothing
    outside (F25's third review: the ) of `case a in a) :;; esac` closed
    the substitution around it), and then reads its text as a command of
    its own, whose here-documents' bodies lie inside the span, as in
    x=`cat <<EOF ... EOF`; so the text is scanned as one, a level deeper
    (scan)."""
    if command[i] != "`":
        return None
    match = BACKQUOTED.match(command, i)
    end, inner = (match.end(), command[i + 1:match.end() - 1]) if match else (len(command), command[i + 1:])
    found = scan(inner, depth=depth + 1)[0] if "<<" in inner else []
    return end, [Body(body.start + i + 1, body.end + i + 1, body.runs) for body in found]


# What closes each opener matching_parenthesis counts: a (, the { of a ${
# and the [ of a $[ or of a [ inside one.
CLOSERS = {"(": ")", "{": "}", "[": "]"}


def opener(ch, dollar, inside):
    """Which opener ch is, where dollar says a $ stands before it, inside
    the opener the count stands in, as parse_matched_pair reads them: a (
    in a parenthesis or after a $, a { after a $ outside $[ ], a [ after a
    $ or inside $[ ]; None for any other character."""
    if ch == "(" and (dollar or inside == "("):
        return "("
    if ch == "{" and dollar and inside != "[":
        return "{"
    if ch == "[" and (dollar or inside == "["):
        return "["
    return None


def matching_parenthesis(command, open_at, memo):
    """Where the ) that closes the ( at command[open_at] stands, counting the
    parentheses after it and stepping over a backslash pair, a quoted span
    and a backquoted one, and reading a ${ to the } that closes it, inside
    which a ( or ) is text but for the ( of a $(, as bash's
    parse_matched_pair does (F25's third review: the ) of ${y:-)} was
    counted), and a $[ to the ] that closes it, a [ inside it counted and
    a ${ not read, $$ opening nothing (its fourth review); None where none
    closes it. A # opens no comment here, since
    parse_matched_pair counts the parentheses after one as it tells
    arithmetic from a subshell: bash reads ((echo a #(, a newline and )))
    as arithmetic. Every ( the count meets outside quotes is matched in the
    same pass and kept in memo, the match the count from it would find, and
    so is every ${ and $[, so a text of many (( that never close is read
    once and not once for each (F25), and so is one where the counts from
    those (( stand in a ${ that never closes (its fourth review: the
    delimiter word X${ of <<X${, which the scan reads as a word and the
    count as a ${, took the count from each (( to the end)."""
    if open_at in memo:
        return memo[open_at]
    opened, i, after_dollar = [("(", open_at)], open_at + 1, False
    others, seen, floor = 0, [], None
    while i < len(command) and opened:
        ch, dollar, after_dollar = command[i], after_dollar, False
        inside = opened[-1][0]
        if ch == "\\":
            i += 2
        elif ch in "'\"":
            i = quote_end(command, i) or len(command)
        elif ch == "`":
            i = backquote_end(command, i)
        elif kind := opener(ch, dollar, inside):
            key = i if kind == "(" else (kind, i)
            if kind != "(" and not others:
                drop = memo.get(("drop", i))
                if drop is not None and len(opened) > drop:
                    floor = len(opened) - drop
                    break
                seen.append((i, len(opened)))
            if key not in memo:
                opened.append((kind, i))
                others += kind != "("
                i += 1
            elif memo[key] is None:
                floor = len(opened)
                break
            else:
                i = memo[key] + 1
        elif ch == CLOSERS[inside]:
            kind, at = opened.pop()
            others -= kind != "("
            memo[at if kind == "(" else (kind, at)] = i
            if not others:
                seen.append((None, len(opened)))
            i += 1
        else:
            # $$, the shell's process id, is a parameter whole, so the {, [
            # or ( after it opens nothing (F25's fourth review).
            after_dollar = ch == "$" and not dollar
            i += 1
    if opened:
        never_closed(memo, opened, seen, len(opened) if floor is None else floor)
    return memo[open_at]


def never_closed(memo, opened, seen, floor):
    """Keeps what a count that ran out of text learned: each opener below
    floor in the stack it left never closes, and from each ${ or $[ it met
    with only parentheses open the count fell no further than to floor or
    to the lowest depth it came back to after it, so a later count that
    meets that opener holding more parentheses than that fall stops there,
    its own ( never closing (F25's fifth review: a text of }<<'X'${y:-((
    units, whose (( stand inside a ${ for the count from each of them,
    took each count to the end)."""
    for kind, at in opened[:floor]:
        memo[at if kind == "(" else (kind, at)] = None
    lowest = floor
    for at, depth in reversed(seen):
        lowest = min(lowest, depth)
        if at is not None:
            memo.setdefault(("drop", at), depth - lowest)


def arithmetic_end(command, i, memo=None):
    """Where the arithmetic whose (( ends just before i ends, past the ))
    that closes it, counting the parentheses inside and stepping over
    quotes and backslash pairs (matching_parenthesis); None where the (( is
    no arithmetic. bash reads (( as arithmetic only where the ) that closes
    its second ( is followed at once by another, and else as a subshell in
    a subshell, as in ((cd x); cat <<EOF ...) or ((echo "))"); ...), whose
    here-documents it reads as any other (F25), and so does this module,
    and where no ) closes it, which bash cannot read, the text after it is
    read as commands. memo keeps the matches across the calls of one scan."""
    close = matching_parenthesis(command, i - 1, {} if memo is None else memo)
    return close + 2 if close is not None and command.startswith(")", close + 1) else None


def delimiter_word(command, i):
    """The word at command[i] that ends a here-document's body, its quotes
    and backslashes taken off as bash takes them off, whether any part of it
    was quoted, and where it ends: at the first character outside quotes
    that ends a word (WORD_BREAKS). None for no word at all, an empty one
    being a word only where quoted. bash reads the word as it reads any, but
    expands nothing in it: a $( ), $(( )), ${ } or backquoted span is part
    of the word as written, a ( or a blank inside it ending nothing, a
    $'...' gives the text it decodes to and a $"..." its text, each quoted
    (F29: bash 5.2 ends the body of <<E$(x)F at the line E$(x)F and that of
    <<$'E\\x4fF' at EOF, where this module ended the word at the ( and read
    the $' as a $)."""
    word, quoted = [], False
    while i < len(command) and command[i] not in WORD_BREAKS:
        ch = command[i]
        if ch == "\\":
            escaped = command[i + 1:i + 2]
            word.append("" if escaped == "\n" else escaped)
            i, quoted = i + 2, quoted or escaped != "\n"
        elif command.startswith("$'", i):
            end = ansi_c_close(command, i + 2)
            word.append(ansi_c_text(command[i + 2:end - 1] if end else command[i + 2:]))
            i, quoted = end or len(command), True
        elif locale_quoted(command, i):
            i += 1
        elif ch in "$`" and (end := expansion_in_word(command, i)):
            word.append(command[i:end])
            i = end
        elif ch in "'\"":
            end = quote_end(command, i)
            inner = command[i + 1:end - 1] if end else command[i + 1:]
            word.append(inner if ch == "'" else DOUBLE_QUOTED_ESCAPE.sub(
                lambda pair: "" if pair.group(1) == "\n" else pair.group(1), inner))
            i, quoted = end or len(command), True
        else:
            word.append(ch)
            i += 1
    text = "".join(word)
    return (text if text or quoted else None), quoted, i


def expansion_in_word(command, i):
    """Where the expansion that opens at command[i] ends, read as bash reads
    one inside a word, a $( ) or $(( )) to the ) that closes it
    (matching_parenthesis), a ${ } to the } that closes it (closing_brace)
    and a backquoted span to its closing ` (backquote_end), each at the end
    of the text where nothing closes it; None where none opens there."""
    if command.startswith("$(", i):
        close = matching_parenthesis(command, i + 1, {})
        return len(command) if close is None else close + 1
    if command.startswith("${", i):
        return closing_brace(command, i + 1)
    if command[i] == "`":
        return backquote_end(command, i)
    return None


def closing_brace(command, i):
    """Where the { at command[i] is closed, past the } that closes it, the
    braces between counted and a quoted span or a backslash pair stepped
    over, as bash's parse_matched_pair reads a ${ }; the end of the text
    where none closes it."""
    depth = 0
    while i < len(command):
        ch = command[i]
        if ch == "\\":
            i += 2
            continue
        if ch in "'\"":
            i = quote_end(command, i) or len(command)
            continue
        depth += (ch == "{") - (ch == "}")
        i += 1
        if not depth:
            return i
    return len(command)


def here_operator(command, i, contexts):
    """Where the scan goes on past the << at command[i], and the
    here-document it opens: past a <<< here-string, which opens none, past
    a << the shell reads as text (takes_here_document), and past the word
    after a << or <<- and the blanks before it, which is the
    here-document's delimiter; None where it opens none, as where no word
    follows, which is a syntax error."""
    if command.startswith("<<<", i):
        return i + 3, None
    if not takes_here_document(contexts):
        return i + 2, None
    strip_tabs = command.startswith("<<-", i)
    start = i + 2 + strip_tabs
    while start < len(command) and command[start] in " \t":
        start += 1
    word, quoted, end = delimiter_word(command, start)
    return end, None if word is None else HereDocument(word, strip_tabs, not quoted, level(contexts))


class Step(NamedTuple):
    """Where the scan goes on after one step, the contexts it then stands
    in, whether a word starts there (a # there opens a comment) and whether
    that word is a command's first (a case there is a case command)."""

    at: int
    contexts: object
    word_start: bool
    command_start: bool


def arithmetic_opened(at, contexts, end):
    """The step into arithmetic whose text starts at command[at] and closes
    at end: a context of its own (ARITHMETIC), read to its closing
    parentheses, since bash runs the substitutions inside it (F29)."""
    return Step(at, Context(ARITHMETIC, contexts, at, end), False, False)


def expansion_step(command, i, contexts, memo):
    """The step of the scan from command[i] where an expansion starts there,
    in double quotes, a parameter expansion or neither, else None: $$, the
    shell's process id, stepped over whole, so the { or [ after it opens
    nothing (F25's fourth review: $${ opened a parameter expansion bash
    never opens, which hid the << after it); an arithmetic expansion opened
    as a context of its own (ARITHMETIC), read to its closing parentheses,
    since bash runs the substitutions inside it, and the here-documents of
    those are told too (F29: until then it was jumped over, and the body of
    a here-document in a substitution in arithmetic was read as commands, a
    quote in it hiding the push after it); or a parameter expansion, old-style arithmetic, $[ ], or a
    substitution opened (F25's third review: a ) inside ${y:-)} closed the
    substitution around it; its fourth: and inside $[ ) ]). The scan reads
    a backquoted span before any step (backquoted)."""
    if command.startswith("$$", i):
        return Step(i + 2, contexts, False, False)
    end = arithmetic_end(command, i + 3, memo) if command.startswith("$((", i) else None
    if end:
        return arithmetic_opened(i + 3, contexts, end)
    if command.startswith("${", i) or command.startswith("$[", i):
        return Step(i + 2, Context(command[i:i + 2], contexts, i + 2), False, False)
    if command.startswith("$(", i):
        return Step(i + 2, Context("$(", contexts), True, True)
    return None


def double_quoted_step(command, i, contexts, memo):
    """One step of the scan inside double quotes, from command[i]: a
    backslash pair, the closing quote, an expansion (expansion_step), or
    one character."""
    if command[i] == "\\":
        return Step(i + 2, contexts, False, False)
    if command[i] == '"':
        return Step(i + 1, contexts.outer, False, False)
    return expansion_step(command, i, contexts, memo) or Step(i + 1, contexts, False, False)


def brace_state(command, i, contexts):
    """The state of the parameter expansion contexts, a ${ }, at command[i],
    as bash's parse_matched_pair keeps it: "first" before its first
    character, "param" in its name, "op" in its operator, "quote" after an
    operator that removes or replaces a pattern, # or % or / or ^ or , after
    the name, and "word" after any other; read from where it was last read
    up to i, and no further once it is settled, so each character is read
    for it once."""
    while contexts.seen < i and contexts.state in ("first", "param", "op"):
        ch, state = command[contexts.seen], contexts.state
        if state == "param" and ch in PATTERN_OPERATORS:
            contexts.state = "quote"
        elif state in ("first", "param") and ch in BRACE_OPERATORS:
            contexts.state = "op"
        elif state == "first":
            contexts.state = "param"
        elif state == "op" and ch not in BRACE_OPERATORS:
            contexts.state = "word"
        contexts.seen += 1
    if contexts.state in ("first", "param", "op"):
        contexts.seen = max(contexts.seen, i)
    return contexts.state


def braced_step(command, i, contexts, memo):
    """One step of the scan inside a parameter expansion, ${ }, from
    command[i]: a backslash pair; a single quote, which outside double
    quotes opens a quoted span and, in double quotes, one after an operator
    that removes or replaces a pattern (brace_state), while after any other
    it opens a span that pairs with the next but in which bash expands
    $( ), backquotes and arithmetic (APOSTROPHE; F24's first review); a
    $'...' outside double quotes, stepped over whole (ansi_c_end; F24's
    first review: ${v:-$'a\\'b'} hid a backquoted push after it); a double
    quote opened, an expansion (expansion_step), the } that closes it, or
    one character, a ( or ) among them, as bash's parse_matched_pair reads
    it."""
    ch = command[i]
    if ch == "\\":
        return Step(i + 2, contexts, False, False)
    if ch == "'":
        if contexts.quoted and brace_state(command, i, contexts) != "quote":
            return Step(i + 1, Context(APOSTROPHE, contexts), False, False)
        return Step(quote_end(command, i) or len(command), contexts, False, False)
    if ch == "$" and not contexts.quoted and command.startswith("$'", i):
        return Step(ansi_c_end(command, i + 2), contexts, False, False)
    if ch == '"':
        return Step(i + 1, Context('"', contexts), False, False)
    if ch == "}":
        return Step(i + 1, contexts.outer, False, False)
    return expansion_step(command, i, contexts, memo) or Step(i + 1, contexts, False, False)


def bracketed_step(command, i, contexts, memo):
    """One step of the scan inside old-style arithmetic, $[ ], or a [ inside
    it, from command[i], as bash's parse_matched_pair reads it: a backslash
    pair, a single quote, which pairs with the next, so a ] between them
    closes nothing, but quotes nothing, since bash expands arithmetic's
    text as it does double quotes' and runs the push in
    echo $[ '$(git push ...)' ] (APOSTROPHE; F24's first review), a
    double quote opened, a [ that the ] after it closes, the ] that closes
    it, a substitution or arithmetic as elsewhere (expansion_step) but for
    a ${, which bash does not read there, so $[ 1 ${x:-]} ] ends at the ]
    inside the braces; or one character, a ( or ) among them."""
    ch = command[i]
    if ch == "\\":
        return Step(i + 2, contexts, False, False)
    if ch == "'":
        return Step(i + 1, Context(APOSTROPHE, contexts), False, False)
    if ch == '"':
        return Step(i + 1, Context('"', contexts), False, False)
    if ch == "[":
        return Step(i + 1, Context("[", contexts), False, False)
    if ch == "]":
        return Step(i + 1, contexts.outer, False, False)
    step = None if command.startswith("${", i) else expansion_step(command, i, contexts, memo)
    return step or Step(i + 1, contexts, False, False)


def apostrophe_step(command, i, contexts, memo):
    """One step of the scan inside a single quote that pairs but quotes
    nothing (APOSTROPHE), from command[i]: the ' that closes it, an
    expansion (expansion_step), or one character, a backslash among them,
    which escapes nothing there (bash 5.2 ends "${v:-'\\'...'}" at the
    second apostrophe)."""
    if command[i] == "'":
        return Step(i + 1, contexts.outer, False, False)
    return expansion_step(command, i, contexts, memo) or Step(i + 1, contexts, False, False)


def arithmetic_step(command, i, contexts, memo):
    """One step of the scan inside arithmetic read apart (ARITHMETIC), from
    command[i]: the )) that closes it, which the count of its parentheses
    found (arithmetic_end, contexts.end), a backslash pair, a single quote,
    which pairs with the next but quotes nothing, since bash expands
    arithmetic's text as it does double quotes' and runs the push in
    echo $(( '$(git push ...)' )) (APOSTROPHE; F24's first review), a
    double quote opened, an expansion (expansion_step), or one character."""
    ch = command[i]
    if i >= contexts.end - 2:
        return Step(contexts.end, contexts.outer, False, False)
    if ch == "\\":
        return Step(i + 2, contexts, False, False)
    if ch == "'":
        return Step(i + 1, Context(APOSTROPHE, contexts), False, False)
    if ch == '"':
        return Step(i + 1, Context('"', contexts), False, False)
    return expansion_step(command, i, contexts, memo) or Step(i + 1, contexts, False, False)


def is_word(command, i, word):
    """Whether the word at command[i] is word, whole."""
    end = i + len(word)
    return command.startswith(word, i) and (end == len(command) or command[end] in WORD_BREAKS)


def case_step(command, i, contexts, word_start, command_start):
    """One step of the scan through a case command inside a substitution,
    or None where the step is any other: bash reads a pattern's ) as
    closing no substitution (F25: x="$(case a in a) cat <<'EOF' ...)" told
    the body at the wrong level). A case among a command's first words
    opens one; its in starts the patterns, where a ( before a pattern opens
    nothing, a ) ends the pattern and esac closes the case; and in a
    pattern's commands ;;, ;& or ;;& starts the next pattern and an esac
    among a command's first words closes it. A ) there, which bash
    refuses, closes the case and what holds it."""
    kind, ch = contexts.kind, command[i]
    if kind == "case" and word_start and is_word(command, i, "in"):
        contexts.kind = "in"
        return Step(i + 2, contexts, False, False)
    if kind == "in":
        if word_start and is_word(command, i, "esac"):
            return Step(i + 4, contexts.outer, False, False)
        if ch in "()":
            contexts.kind = "case)" if ch == ")" else kind
            return Step(i + 1, contexts, True, ch == ")")
        return None
    if kind == "case)":
        if command.startswith(";;", i) or command.startswith(";&", i):
            contexts.kind = "in"
            return Step(i + (3 if command.startswith(";;&", i) else 2), contexts, True, False)
        if word_start and command_start and is_word(command, i, "esac"):
            return Step(i + 4, contexts.outer, False, False)
        if ch == ")":
            return Step(i + 1, contexts.outer.outer, True, False)
    if word_start and command_start and is_word(command, i, "case"):
        return Step(i + 4, Context("case", contexts), False, False)
    return None


def unquoted_step(command, i, contexts, word_start, command_start, memo):
    """One step of the scan outside double quotes, at the top or inside a
    substitution, from command[i] (Step). A backslash pair, a single-quoted
    span and a comment are jumped over, a comment to its newline; a double
    quote, arithmetic, (( )) or $(( )), a parameter expansion and a
    substitution open a context (expansion_step), a parenthesis inside a
    substitution one the ) after it closes, and a case inside one is read
    as bash reads it (case_step). A line continuation leaves whether a word
    and a command start where they were, since bash takes it out before it
    reads on (F29: a # after one that follows a blank opens a comment)."""
    ch = command[i]
    if ch == "\\":
        joins = joins_lines(command[i:i + 2])
        return Step(i + 2, contexts, word_start and joins, command_start and joins)
    if ch == "'":
        return Step(quote_end(command, i) or len(command), contexts, False, False)
    if ch == "$" and command.startswith("$'", i):
        return Step(ansi_c_end(command, i + 2), contexts, False, False)
    if ch == "#" and word_start:
        end = command.find("\n", i)
        return Step(len(command) if end == -1 else end, contexts, False, command_start)
    end = arithmetic_end(command, i + 2, memo) if command.startswith("((", i) else None
    if end:
        return arithmetic_opened(i + 2, contexts, end)
    if ch == '"':
        return Step(i + 1, Context('"', contexts), False, False)
    step = expansion_step(command, i, contexts, memo)
    if step:
        return step
    step = case_step(command, i, contexts, word_start, command_start) if contexts is not None else None
    if step:
        return step
    if ch == "(" and contexts is not None:
        return Step(i + 1, Context("(", contexts), True, True)
    if ch == ")" and contexts is not None:
        return Step(i + 1, contexts.outer, True, False)
    if word_start and command_start and ch in BEFORE_A_COMMAND_OPENS:
        word = next((word for word in BEFORE_A_COMMAND if is_word(command, i, word)), None)
        if word:
            return Step(i + len(word), contexts, False, True)
    return Step(i + 1, contexts, ch in WORD_BREAKS, ch in ";&|(\n" or (command_start and ch in " \t"))


def scan_step(command, i, contexts, word_start, command_start, memo):
    """One step of the scan from command[i], inside double quotes, a
    parameter expansion, old-style arithmetic, a single quote that quotes
    nothing, arithmetic read apart or none of them (double_quoted_step,
    braced_step, bracketed_step, apostrophe_step, arithmetic_step,
    unquoted_step)."""
    if contexts is not None and contexts.kind == '"':
        return double_quoted_step(command, i, contexts, memo)
    if contexts is not None and contexts.kind == "${":
        return braced_step(command, i, contexts, memo)
    if contexts is not None and contexts.kind in ("$[", "["):
        return bracketed_step(command, i, contexts, memo)
    if contexts is not None and contexts.kind == APOSTROPHE:
        return apostrophe_step(command, i, contexts, memo)
    if contexts is not None and contexts.kind == ARITHMETIC:
        return arithmetic_step(command, i, contexts, memo)
    return unquoted_step(command, i, contexts, word_start, command_start, memo)


def leave(pending, before, after):
    """The documents still pending at each level the scan left in going from
    the contexts before to those after, moved to the level after stands at,
    after those already pending there: bash reads the body of a
    here-document whose substitution or double quotes close on its
    operator's own line from the next newline it reads outside them (F25:
    x=$(cat <<'EOF') read its body from the line after the ), and the scan
    waited for a newline inside the substitution, which never came)."""
    if after is not None and after.outer is before:
        return
    context, target = before, level(after)
    while context is not None and context is not after:
        if context.level is context and context in pending:
            moved = [document._replace(within=target) for document in pending.pop(context)]
            pending.setdefault(target, []).extend(moved)
        context = context.outer


def scan(command, i=0, contexts=None, depth=0, memo=None):
    """The bodies of the here-documents the scan meets in command from i,
    standing in contexts, in order, and where it stops: past the ) that
    closes the substitution it starts in, or at the end of the text. A body
    starts at the first newline after its operator that the shell reads at
    the operator's own level (level), past a quote that spans lines, or at
    the level it is moved to where that level closes first (leave), and
    ends past the line that is its delimiter (body_end). A << is an
    operator outside quotes and in a substitution, double-quoted or not
    (takes_here_document), and none in a <<< here-string, in arithmetic, in
    a comment or in a quoted span, which the scan steps over as the shell
    reads them (scan_step), and a backquoted span is read whole, its own
    here-documents told inside it (backquoted). Each step takes the same
    time however deep the scan stands, and each ( is matched once
    (matching_parenthesis) through
    memo, kept for the text across every scan of it, so the scan takes a
    time in proportion to the text, times the depth its bodies nest at
    (depth, MAX_NESTING). A memo kept per scan, as until F25's third
    review, matched a ( again in the scan of each substitution of an
    expanded body, a time in the square of the body's length. The reading
    stops where it runs past its deadline (in_time)."""
    if depth > MAX_NESTING:
        raise TooDeep()
    closes = contexts is not None
    bodies, pending, word_start, command_start = [], {}, True, True
    memo = {} if memo is None else memo
    steps = 0
    while i < len(command):
        steps += 1
        if not steps % CHECKS_EVERY:
            in_time()
        plain = plain_run(command, i, contexts, word_start)
        if plain:
            i, word_start, command_start = plain, False, False
            continue
        if command[i] == "\n" and level(contexts) in pending:
            found = bodies_from(command, i + 1, pending.pop(level(contexts)), depth)
            bodies += found
            i, word_start, command_start = found[-1].end, True, True
        elif command.startswith("<<", i):
            i, document = here_operator(command, i, contexts)
            if document:
                pending.setdefault(document.within, []).append(document)
            word_start = command_start = False
        elif (span := backquoted(command, i, depth)) is not None:
            i, inner = span
            bodies += inner
            word_start = command_start = False
        else:
            before = contexts
            i, contexts, word_start, command_start = scan_step(command, i, contexts, word_start, command_start, memo)
            if contexts is not before:
                leave(pending, before, contexts)
                if closes and contexts is None:
                    return bodies, i
    return bodies, len(command)


def substitutions_in(command, start, end, depth=0):
    """The command substitutions, each $( ) whole, in the text of a body the
    shell expands, from start to end: the shell runs each when it expands
    the body, and none of the body's other text, whose quotes are text too
    (F25: "It's $(git push ...) done" hid the push between two
    apostrophes). bash reads the body's text alone, so each substitution is
    read inside it: it ends where the scan of it closes it (scan), and what
    is kept of it is its text with its own here-documents' bodies cut out
    in turn (cut_out), since their quotes are text as well (F25: a body in
    a substitution in a body, whose apostrophe paired with one after the
    push). A backquoted span is a substitution too, read with its own
    bodies cut out (backquoted; F25's fourth review). A backslash keeps the
    $ after it, and $(( )) is arithmetic, kept whole where a substitution
    stands inside it, since bash runs it, which the reading of the command's
    substitutions then finds (substituted; F24's first review: a body line
    echo "$(( $(git push ...) ))" ran the push). One
    memo serves every scan of the text, so each ( in it is matched once
    (scan)."""
    text, found, memo, i = command[start:end], [], {}, 0
    while (mark := BODY_MARK.search(text, i)) is not None:
        in_time()
        i = mark.start()
        arithmetic = arithmetic_end(text, i + 3, memo) if text.startswith("$((", i) else None
        if text[i] == "\\":
            i += 2
        elif arithmetic:
            if "$(" in text[i + 3:arithmetic] or "`" in text[i + 3:arithmetic]:
                found.append(text[i:arithmetic])
            i = arithmetic
        elif text.startswith("$(", i):
            bodies, close = scan(text, i + 2, Context("$("), depth + 1, memo)
            found.append(cut_out(text, bodies, i, close))
            i = close
        elif (span := backquoted(text, i, depth)) is not None:
            found.append(cut_out(text, span[1], i, span[0]))
            i = span[0]
        else:
            i += 1
    return tuple(found)


# What substitutions_in stops at in an expanded body: a backslash, a $ and
# a backquote, the only characters that open or hide a substitution there,
# so the rest of a long body is passed over at once (F24's first review).
BODY_MARK = re.compile(r"[\\$`]")


def closes_substitution(line, delimiter):
    """Whether a line of a body inside a substitution ends the body though
    it is not the delimiter: one that starts with the delimiter and holds a
    ) after it, as "EOF)" does, where bash 5.2 ends the body after the
    delimiter, with a warning, and reads the rest of the line as the
    substitution's own text."""
    return line.startswith(delimiter) and ")" in line[len(delimiter):]


def continues(text):
    """Whether a line of a body the shell expands goes on into the next: it
    ends in a backslash that no backslash before it escapes, an odd run of
    them, which bash takes out with the newline after it before it compares
    the line with the delimiter (F25: a line of E and a backslash, then one
    of OF, ends a body whose delimiter is EOF)."""
    return (len(text) - len(text.rstrip("\\"))) % 2 == 1


def body_line(command, i, joins):
    """The line of a body that starts at command[i] as bash compares it with
    the delimiter, the spans of command it is made of, and where the
    newline that ends it stands: in a body the shell expands (joins), a
    line goes on into the next while it continues (continues), and in a
    quoted one each line stands alone."""
    pieces, spans = [], []
    while True:
        newline = command.find("\n", i)
        line_end = len(command) if newline == -1 else newline
        if not (joins and newline != -1 and continues(command[i:line_end])):
            pieces.append(command[i:line_end])
            spans.append((i, line_end))
            return "".join(pieces), spans, line_end
        pieces.append(command[i:line_end - 1])
        spans.append((i, line_end - 1))
        i = line_end + 1


def place_in(spans, k):
    """Where in the command the character k of a line made of spans
    (body_line) stands."""
    for start, end in spans:
        if k < end - start:
            return start + k
        k -= end - start
    return spans[-1][1]


def body_end(command, start, document):
    """Where a here-document's body that starts at start ends: past the
    first line that is its delimiter (body_line), that line's leading tabs
    taken off after <<-, or, inside a substitution, right after the
    delimiter on a line that ends the body there (closes_substitution); and
    at the end of the text, where bash ends a body no line ends."""
    i, inside = start, document.within is not None and document.within.in_substitution
    while i < len(command):
        line, spans, line_end = body_line(command, i, document.expands)
        indent = len(line) - len(line.lstrip("\t")) if document.strip_tabs else 0
        if line[indent:] == document.delimiter:
            return min(line_end + 1, len(command))
        if inside and closes_substitution(line[indent:], document.delimiter):
            return place_in(spans, indent + len(document.delimiter))
        i = line_end + 1
    return len(command)


def bodies_from(command, start, documents, depth=0):
    """The bodies of documents, here-documents whose operators stand on one
    line, read in order from start, the character after the newline that
    ends that line: each starts where the one before it ends (body_end), and
    what of it the shell runs is read where it expands the body
    (substitutions_in)."""
    bodies = []
    for document in documents:
        end = body_end(command, start, document)
        runs = substitutions_in(command, start, end, depth) if document.expands else ()
        bodies.append(Body(start, end, runs))
        start = end
    return bodies


def here_document_bodies(command):
    """The bodies of a shell command's here-documents (Body), in order,
    which the shell reads as text and not as commands (F25), as the scan
    from the top tells them (scan); TooDeep where they nest deeper than
    MAX_NESTING."""
    return scan(command)[0] if "<<" in command else []


def cut_out(command, bodies, start=0, end=None, runs=True):
    """The command from start to end, by default the whole, with each body,
    in order and apart, taken out, and in its place what of it the shell
    runs, each on a line of its own, or nothing where not runs."""
    kept, at = [], start
    for body in bodies:
        kept += [command[at:body.start], *(run + "\n" for run in body.runs if runs)]
        at = body.end
    return "".join(kept) + command[at:len(command) if end is None else end]


# --- substitutions ---------------------------------------------------------------

# What may open a substitution, so a command that holds none of them is
# not scanned for one.
SUBSTITUTION_MARKS = ("$(", "`", "<(", ">(", "((")
# The characters a substitution other than a backquoted span opens at
# (substitution_opens), so the scan asks only there.
OPENERS = set("$<>(")
# A run of characters no step of the scan reads but as part of a word,
# inside a word: the scan of substitutions steps over it at once, since a
# step there only moves on (scan_step), and most of a long command is such
# runs.
PLAIN_RUN = re.compile(r"[^\s$`'\"\\#(){}\[\]<>;&|]+")
# The words a step of the scan reads as more than part of a word where one
# starts (unquoted_step, case_step), so a plain run that is one of them is
# read a character at a time (plain_run).
STEPPED_WORDS = frozenset(BEFORE_A_COMMAND + ("case", "in", "esac"))
# The contexts in which a word's start means more than in any other: a case
# command's, whose patterns and words case_step reads.
CASE_KINDS = ("case", "in", "case)")


def plain_run(command, i, contexts, word_start):
    """Where the run of plain characters (PLAIN_RUN) at command[i] ends, or
    None where none starts there or one must be read a character at a time:
    at a word's start, a run that is one of STEPPED_WORDS or stands in a
    case command (CASE_KINDS). A step through such a run only moves on, so
    both scans step over it at once (F24's first review: a body of 52,000
    $(a) took 3.7 seconds read a character at a time)."""
    run = PLAIN_RUN.match(command, i)
    if run is None or not word_start:
        return run and run.end()
    if (contexts is not None and contexts.kind in CASE_KINDS) or run.group(0) in STEPPED_WORDS:
        return None
    return run.end()


# A run of blanks or newlines outside quotes and expansions, where the scan
# of substitutions only starts a word, and after a newline a command too
# (blank_run).
BLANKS = re.compile(r"[ \t]+")
NEWLINES = re.compile(r"[ \t\n]*\n[ \t]*")
# The contexts in which a blank or a newline only starts a word or a
# command: the top, a substitution and a parenthesis inside one.
READS_BLANKS = (None, "$(", "(")


def blank_run(command, i, contexts, command_start):
    """Where the run of blanks at command[i] ends and whether a command
    starts after it, or None where none starts there or contexts is one in
    which a blank means more (READS_BLANKS): blanks keep whether a command
    starts, and a newline among them starts one, as unquoted_step reads
    each."""
    if (contexts.kind if contexts is not None else None) not in READS_BLANKS:
        return None
    run = NEWLINES.match(command, i)
    if run:
        return run.end(), True
    run = BLANKS.match(command, i)
    return (run.end(), command_start) if run else None


# How many steps the scan of substitutions takes between two looks at the
# deadline (in_time): a look costs about as much as a step, and a step
# takes a few microseconds, so the deadline is still kept to well under a
# millisecond.
CHECKS_EVERY = 64
# A backslash pair inside a backquoted span that bash takes the backslash
# off before it runs the span's text: before $, ` and a backslash, and,
# where the span stands in double quotes, before a double quote too (bash
# 5.2 runs git push origin main for "`git push origin \"main\"`").
BACKQUOTE_ESCAPE = re.compile(r"\\([$`\\])")
DOUBLE_QUOTED_BACKQUOTE_ESCAPE = re.compile(r'\\([$`\\"])')


class Substitution:
    """A substitution's own text as it is read apart (substitutions): the
    pieces of it read so far, each nested substitution standing as
    SUBSTITUTED; where in the scanned text the next piece starts, None
    while a substitution nested in it is open, which holds the rest of the
    text where it never closes, so no text is read in two of them; the
    substitutions nested in it, in order; and the context the scan opened
    for it, whose closing ends it (None for a backquoted span or the text
    itself, read whole)."""

    __slots__ = ("pieces", "at", "nested", "context")

    def __init__(self, at=0, context=None):
        self.pieces, self.at, self.nested, self.context = [], at, [], context


def reads_commands(contexts):
    """Whether the scan of substitutions stands where bash reads commands
    in contexts: outside double quotes, a parameter expansion and
    old-style arithmetic, at the top or in a substitution or a subshell,
    where a <( or >( opens a process substitution and a (( arithmetic; in
    double quotes each is text."""
    return contexts is None or contexts.kind not in TEXT_ONLY


def substitution_opens(command, i, contexts, memo):
    """Where the own text of the substitution that opens at command[i]
    starts, and, for arithmetic, where it ends, past its )) (arithmetic_end),
    or None where none opens there: a $( anywhere the scan reads, in double
    quotes, a parameter expansion and arithmetic too, and $(( )) among
    them, since bash runs the substitutions inside arithmetic; and, where
    bash reads commands (reads_commands), a <( or >( and arithmetic,
    (( )). Arithmetic is read in a context of its own (arithmetic_step),
    where a single quote quotes nothing, so the substitutions inside it
    are found, which bash runs, and its own text, read as commands, runs
    nothing but a word with no push in it."""
    if command.startswith("$((", i) and (end := arithmetic_end(command, i + 3, memo)):
        return i + 3, end
    if command.startswith("$(", i):
        return i + 2, None
    if not reads_commands(contexts):
        return None
    if command.startswith("<(", i) or command.startswith(">(", i):
        return i + 2, None
    if command.startswith("((", i) and (end := arithmetic_end(command, i + 2, memo)):
        return i + 2, end
    return None


def in_double_quotes(contexts):
    """Whether the scan stands directly in double quotes at contexts."""
    return contexts is not None and contexts.kind == '"'


def closed_by(contexts, substitution):
    """Whether the contexts a step of the scan left it in are outside the
    one opened for substitution, so the step read the ) that closes it: a
    step closes one context, or a case pattern's ) the case and the
    substitution that holds it (case_step), so the contexts are the
    substitution's own outer one."""
    opened = substitution.context
    return opened is not None and contexts is opened.outer


def backquoted_substitution(command, i, quoted=False):
    """The backquoted span that opens at command[i] read apart, as bash
    runs its text, a backslash before $, ` or a backslash taken off, and
    before a double quote where the span stands in double quotes (quoted),
    the bodies of the here-documents that text then holds cut out
    (cut_out) where the span as written holds none, which the reading cut
    out already (readings), since bash reads them as text (F24's first
    review's fuzz:
    `"\\$(cat <<'E' ...)"` hid a here-document from the scan of the
    command, whose body's quote then paired with the push's line), its
    own substitutions read inside it in turn (substitutions), and
    where the scan goes on past it; None where no span opens there. The
    span runs to the end of the text where no ` closes it. bash nests a
    backquoted span in another only by backslashes, which double at each
    level, so the depth stays below the logarithm of the length."""
    if command[i] != "`":
        return None
    end = backquote_end(command, i)
    inner = command[i + 1:end - 1] if BACKQUOTED.match(command, i) else command[i + 1:]
    escape = DOUBLE_QUOTED_BACKQUOTE_ESCAPE if quoted else BACKQUOTE_ESCAPE
    text = escape.sub(r"\1", inner)
    if "<<" in text and not here_document_bodies(inner):
        text = cut_out(text, here_document_bodies(text))
    return substitutions(text), end


def ansi_c_end(command, i):
    """Where bash's $'...' that opened before command[i] ends (ansi_c_close),
    or at the end of the command where none closes it. Read as an ordinary
    quoted span, $'it\\'s' closed at its escaped apostrophe, and the
    apostrophe after it opened a quote that hid every substitution up to the
    next one (F24's fuzz)."""
    return ansi_c_close(command, i) or len(command)


def substitutions(command):
    """The command as a Substitution read whole, its own text the command
    with each substitution in it standing as SUBSTITUTED and each read
    apart, nested in it as the shell nests them (F24): a $( ), a backquoted
    span, a <( ) or >( ) and arithmetic, as substitution_opens tells them,
    in double quotes and a parameter expansion too, but never in single
    quotes, a comment or a backslash pair, which the scan steps over as the
    shell reads them (scan_step), nor in a $'...' where the shell reads
    commands or in a parameter expansion outside double quotes, stepped
    over whole with its escapes (ansi_c_end), nor in a here-document's
    delimiter word, which bash reads as text (here_operator; F29: the
    $(echo x) of <<E$(echo x)F was read as a command bash runs), each
    running to the ) that closes it as
    the scan of here-documents reads one (closed_by). The command is read
    with its here-documents' bodies cut out (readings), so a quoted body,
    as a commit message's, is text here too. The scan is one pass, each
    step in the same time however deep it stands, so the reading takes a
    time in proportion to the text, and a text that holds no
    SUBSTITUTION_MARKS is its own text whole, unscanned."""
    whole = Substitution()
    if not any(mark in command for mark in SUBSTITUTION_MARKS):
        whole.pieces.append(command)
        return whole
    open_ones, contexts, word_start, command_start, memo, i, steps = [whole], None, True, True, {}, 0, 0
    while i < len(command):
        steps += 1
        if not steps % CHECKS_EVERY:
            in_time()
        plain = plain_run(command, i, contexts, word_start)
        if plain:
            i, word_start, command_start = plain, False, False
            continue
        blanks = blank_run(command, i, contexts, command_start)
        if blanks:
            (i, command_start), word_start = blanks, True
            continue
        current, ch = open_ones[-1], command[i]
        if ch == "<" and command.startswith("<<", i) and reads_commands(contexts):
            i, word_start, command_start = here_operator(command, i, contexts)[0], False, False
            continue
        opens = substitution_opens(command, i, contexts, memo) if ch in OPENERS else None
        if opens is not None:
            start, end = opens
            current.pieces += [command[current.at:i], SUBSTITUTED]
            current.at = None
            contexts = Context("$(", contexts) if end is None else Context(ARITHMETIC, contexts, start, end)
            nested = Substitution(start, contexts)
            current.nested.append(nested)
            open_ones.append(nested)
            i, word_start, command_start = start, True, True
        elif ch == "`" and (span := backquoted_substitution(command, i, in_double_quotes(contexts))) is not None:
            nested, end = span
            current.pieces += [command[current.at:i], SUBSTITUTED]
            current.nested.append(nested)
            current.at, i, word_start, command_start = end, end, False, False
        else:
            i, contexts, word_start, command_start = scan_step(command, i, contexts, word_start, command_start, memo)
            while closed_by(contexts, open_ones[-1]):
                closed = open_ones.pop()
                closed.pieces.append(command[closed.at:i - (2 if closed.context.kind == ARITHMETIC else 1)])
                open_ones[-1].at = i
    for still_open in open_ones:
        if still_open.at is not None:
            still_open.pieces.append(command[still_open.at:])
    return whole


def substituted(command):
    """What the shell runs of the substitutions in a command, as text to
    read after it: each substitution's own text, made plain leniently
    (plain_text), as the commands of a subshell, with those nested in it
    after its own text inside the same subshell, so a cd before them in it
    holds for them, as it does in bash (substitutions); one whose own text
    names no command, nothing in it but the substitutions nested in it and
    characters that make no word (WORD_CHARACTER), runs no command and
    moves no shell, so those are read in its place, which keeps a run of
    nested ones short; and of the substitutions of one text each that
    reads like one before it is read once (alike, distinct). Empty where
    the command holds none. Until F24's fix
    a substitution in double quotes was read as part of its word, so
    echo "$(git push ...)" and out="$(git push ...)" passed, and a <( ) or
    >( ) as a redirection that took the word after it, so
    cat <(git push ...) passed too."""
    return read_apart(command)[1]


def read_apart(command):
    """The command's own text, each substitution in it standing as
    SUBSTITUTED, and what the shell runs of its substitutions
    (substituted); the command itself and nothing where it holds none."""
    if not any(mark in command for mark in SUBSTITUTION_MARKS):
        return command, ""
    root = substitutions(command)
    keys = alike(root)
    out, pending = [], [(nested, False) for nested in reversed(distinct(root.nested, keys))]
    while pending:
        in_time()
        item, done = pending.pop()
        if done:
            out.append(" ) ")
            continue
        own = "".join(item.pieces)
        if WORD_CHARACTER.search(own.replace(SUBSTITUTED, "")):
            out += [" ; ( ", plain_text(own, lenient=True, spans_run=False)]
            pending.append((item, True))
        pending += [(nested, False) for nested in reversed(distinct(item.nested, keys))]
    return "".join(root.pieces), "".join(out)


# What a substitution's own text holds where it may name a command: a
# letter, a digit or an underscore, outside the substitutions nested in it.
# One with none, as "( or a lone quote, names no git, no cd and nothing
# else this module reads, so substituted reads the substitutions nested in
# it in its place (F24's first review: "$( repeated 87,000 times was read
# as as many subshells).
WORD_CHARACTER = re.compile(r"\w")


def alike(root):
    """A key for each substitution under root (Substitution), by identity,
    the same for two that read alike: the same own text, and the same
    substitutions nested in it, told by their keys, each once. Built from
    the innermost out, one step per substitution, so it takes a time in
    proportion to the text however deep they nest."""
    keys, interned, stack = {}, {}, [(root, False)]
    while stack:
        item, ready = stack.pop()
        if not ready:
            stack.append((item, True))
            stack += [(nested, False) for nested in item.nested]
            continue
        inner = tuple(dict.fromkeys(keys[id(nested)] for nested in item.nested))
        keys[id(item)] = interned.setdefault(("".join(item.pieces), inner), len(interned))
    return keys


def distinct(substitutions_here, keys):
    """The substitutions of one text that read unlike any before them
    (alike), in order. Each is read where the shell stands after the text
    they stand in, in a subshell of its own, so two that read alike are
    judged alike, and substituted reads such a text once (F24's first
    review: a body of 52,000 $(a) lexed each)."""
    seen, kept = set(), []
    for nested in substitutions_here:
        if keys[id(nested)] not in seen:
            seen.add(keys[id(nested)])
            kept.append(nested)
    return kept


def readings(command):
    """The texts this module reads a shell command as, each made plain with
    a quote that never closes read leniently (plain_text): the command
    whole, its here-documents' bodies read as commands, as before F25's
    fix; and, where it holds bodies this module tells
    (here_document_bodies), the command with them cut out, since a body's
    quotes are none the shell reads, and what the shell runs of them read
    apart after it with the command's own substitutions (substituted), as
    since F24's first review it is read only there, where it stood in the
    body's place as well until then.
    Every push either finds is judged, so a body misread hides no command
    the whole reading shows. TooDeep where the bodies nest deeper than this
    module reads. What the shell runs of the command's substitutions, a
    backquoted span's text among them, is read apart after the reading
    (substituted; F24), and the whole reading reads it only where the
    command holds no body: where it does, a substitution in a body, as a
    commit message's `git push --force` in Markdown, is text bash never
    runs, while the reading with the bodies cut out reads every
    substitution bash runs, those of an expanded body among them (F25's
    fifth review read backquoted spans so; F24's fix reads every
    substitution so, and plain_text runs no span of its own here).

    And where the command holds substitutions, it is read once more as its
    own text, each substitution standing in its word as SUBSTITUTED
    (read_apart), so a word that a substitution builds a part of, as
    git pu$(true)sh, which the whole reading cuts at the substitution's
    parentheses, is read as the one word bash makes of it, with a hole
    where the substitution stands (word_pieces; F29)."""
    bodies = here_document_bodies(command)
    cut = cut_out(command, bodies) if bodies else command
    own, apart = read_apart(cut)
    if not bodies:
        found = [plain_text(command, lenient=True, spans_run=False) + apart]
    else:
        found = [plain_text(command, lenient=True, spans_run=False),
                 plain_text(cut_out(command, bodies, runs=False), lenient=True, spans_run=False) + apart]
    return found + own_text_reading(own, cut) + pattern_reading(own)


def pattern_reading(own):
    """The reading of a command's own text (read_apart) with each extended
    pattern, as @(main), a part of its word (pattern_group), where it holds
    any: bash reads them so with extglob on, a shopt a session may set, and
    fills one from the names of files, so git push origin @(main) pushes
    main where a file is named so (F29's first review). The readings
    without it stand too, as bash reads the text with extglob off."""
    return [plain_text(own, lenient=True, spans_run=False, patterns=True)] if EXTENDED_PATTERN.search(own) else []


EXTENDED_PATTERN = re.compile(r"(?<!\$)[@?*+!]\(")


def own_text_reading(own, cut):
    """The reading of a command's own text, each substitution a hole in its
    word (read_apart), where it holds any: nothing where own is the command
    itself."""
    return [plain_text(own, lenient=True, spans_run=False)] if own != cut else []


def shell_parts(command):
    """A shell command as the shell's parts (parts_of), read whole and
    strictly: a command shlex cannot read whole, as one with a quote that
    never closes, yields nothing, which simple_commands, the template's
    shim's reading, keeps. The hook's own reading reads such a command
    leniently and more than one way (placed_commands, readings)."""
    return parts_of(plain_shell(command))


def parts_of(text):
    """A text made plain (plain_text) as the shell's parts, in order: each
    simple command as the list of its words, and each run of operators
    between them, such as ;, && or ) |, as a string; nothing where shlex
    cannot read it.

    shlex with punctuation_chars keeps the shell's operators as their own
    tokens: ;, &&, | and ( ) end one command, a redirection swallows its
    target and the file descriptor digit before it, and a parenthesis in a
    redirection's run, as in "(...)>log", still closes or opens its
    subshell, after the command it ends. Newlines and comments are read by
    plain_text first, so shlex's own comment handling, which would read a #
    inside a word as a comment and swallow the rest of the text, is off.
    """
    lexer = shlex.shlex(text, posix=True, punctuation_chars=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    tokens = []
    try:
        for token in lexer:
            if not len(tokens) % CHECKS_EVERY:
                in_time()
            tokens.append(token)
    except ValueError:
        return []
    parts, current, eat_next, held, fed = [], Words(), False, "", None
    for token in tokens:
        # An empty word, '' or "", is a word and no operator: read as one, it
        # ended the command, so "git -C '' push origin main" read as no push
        # until F21's fix.
        operator = bool(token) and all(ch in "();<>|&" for ch in token)
        if operator and (">" in token or "<" in token):
            if current and current[-1].isdigit():
                current.pop()
            eat_next = "<<<" if token.endswith("<<<") else True
            held += "".join(ch for ch in token if ch in "()")
        elif operator:
            if current:
                parts.append(current)
            parts.append(held + token)
            current, eat_next, held = Words(), False, ""
        elif eat_next:
            if eat_next == "<<<":
                current.fed = token
            eat_next = False
        else:
            current.append(token)
    if current:
        parts.append(current)
    if held:
        parts.append(held)
    return parts


class Words(list):
    """A simple command's words (parts_of), with the text a <<< here-string
    feeds it on its standard input, or None (fed): a shell reads that text
    as commands, as bash <<<'git push origin main' runs the push (F32)."""

    fed = None


def simple_commands(command):
    """The argument lists of the simple commands inside a shell command,
    read whole and strictly (shell_parts)."""
    return [part for part in shell_parts(command) if isinstance(part, list)]


def program(word):
    """A command's name without its directory: /usr/bin/git is git."""
    return word.rsplit("/", 1)[-1]


def value_is_next(bundle, value_letters):
    """Whether a bundle of short options ends in one that takes a value, so
    that the value is the next word, as in -n alone or -iu; a value
    attached, as in -n5 or -uHOME, is the rest of the bundle."""
    for at, letter in enumerate(bundle):
        if letter in value_letters:
            return at == len(bundle) - 1
    return False


def long_option_of(word, wrapper):
    """The wrapper's long option a word that starts with -- names, as GNU's
    getopt_long resolves it, the word read up to any =: the option of
    exactly that name, else the one option whose name it starts, as --uns
    names --unset; None where it starts the names of more than one option,
    as env's --d does --default-signal's and --debug's, or of none, either
    of which the wrapper refuses, running nothing. glibc also takes a start
    shared by two names of one option, which no wrapper here has."""
    name = word[2:].partition("=")[0]
    if name in wrapper.long_options:
        return name
    found = [option for option in wrapper.long_options if option.startswith(name)]
    return found[0] if len(found) == 1 else None


def long_value_is_next(word, wrapper):
    """Whether a wrapper's long option word takes the next word as its
    value: one that names an option that needs a value (VALUE), with none
    attached by =. One whose value is OPTIONAL takes it only attached, as
    env's --default-signal does, and a FLAG none. A word long_option_of
    cannot resolve is read as taking none, as every long option the guard
    did not know was before F26's fix, and the words after it as the
    wrapper reads them: the wrapper runs nothing then, so a push after it
    refused costs a rephrasing of a command that runs nothing, and a
    wrapper of another version, which may read the word as one option, as
    an env from before --ignore-signal reads --i as --ignore-environment,
    has its push judged. No start that options of one wrapper share here
    is one of an option that takes the next word, so the word is read as
    each of the options it starts would read it."""
    return "=" not in word and wrapper.long_options.get(long_option_of(word, wrapper)) == VALUE


def ends_options(word):
    """Whether a word ends a wrapper's options where it stands, unread as
    one: a lone -, which getopt reads as no option, as it reads every word
    that does not start with - (F28's review: the guard read on past it, so
    "env - --split-string=x git push" was read as splitting x, where env
    takes the - as -i and --split-string=x as an assignment)."""
    return word == "-"


def takes_next_word(word, wrapper):
    """Whether a wrapper's option word takes the next word as its value
    (long_value_is_next, value_is_next)."""
    if word.startswith("--"):
        return long_value_is_next(word, wrapper)
    return value_is_next(word[1:], wrapper.short_values)


def read_options(rest, read, wrapper, splits=None):
    """Move a wrapper's options, each with its value (Wrapper,
    long_value_is_next), and a -- that ends them, from rest, the words
    still to read, the next last, to read, stopping before the first word
    that is no option or a lone - (ends_options); for env (splits), each
    string -S gives is split and its words read in its place (spliced),
    splits the budget it is split by."""
    while rest and rest[-1][0].startswith("-") and not ends_options(rest[-1][0]):
        in_time()
        word, depth = rest.pop()
        if word == "--":
            read.append(word)
            return
        if splits and spliced(rest, read, word, depth, splits):
            continue
        read.append(word)
        if rest and takes_next_word(word, wrapper):
            read.append(rest.pop()[0])


# What GNU env's -S reads as the space between words, and what a backslash
# gives before each character it takes, besides _ and c (split_string).
SPLIT_SPACES = " \t\n\v\f\r"
SPLIT_ESCAPES = {"f": "\f", "n": "\n", "r": "\r", "t": "\t", "v": "\v", '"': '"', "#": "#", "$": "$", "'": "'",
                 "\\": "\\"}
SPLIT_VARIABLE = re.compile(r"\$\{[A-Za-z_][A-Za-z0-9_]*\}")
# The refusal of a command whose env -S strings nest, a string's words
# holding an -S of their own, more than MAX_NESTING deep, or whose strings
# hold more than SPLIT_BUDGET characters in all, each level split again:
# the guard refuses it unread rather than run out of time (F28's review
# and its second).
SPLIT_BUDGET = 2 * MAX_COMMAND
SPLIT_TOO_DEEP = (f"a command whose env -S strings nest more than {MAX_NESTING} deep, or hold more than "
                  f"{SPLIT_BUDGET} characters to split in all, is refused: the guard cannot read it in time; give "
                  "env the words themselves (.claude/rules/board.md)")
# The most distinct ${NAME}s of a command's -S strings that the guard reads
# each way, set and unset, one reading per subset (unset_subsets), and
# then only while the readings split no more than SPLIT_BUDGET characters
# in all; a command with more that may hold a push is refused unread
# (F28's second review and its third).
MOST_NAMES = 6
SPLIT_TOO_MANY_NAMES = (f"a command whose env -S strings hold more than {MOST_NAMES} distinct ${{NAME}}s, or too "
                        "much text to read once for each mix of them, is refused when its words may spell push: "
                        "env drops each name that is not set, and the guard reads every mix only up to "
                        f"{MOST_NAMES}; give env the words themselves (.claude/rules/board.md)")


def split_string(text, unset=False):
    """The words env's -S splits text into, as GNU coreutils' env.c
    (build_argv) splits them: at spaces outside quotes; in single quotes
    every character as written but \\\\ and \\', in double quotes and
    outside them the backslash's escapes; \\_ a space between words
    outside double quotes and a space in the word inside them, \\c the end
    of the string, and a # that starts a word a comment to the end. A
    ${NAME} is kept as written, as the guard keeps the shell's own $NAME,
    since its value is the one in the command's environment; or, unset
    (True, or a set holding its name), expands to nothing, as env expands a variable that is not set, starting
    no word, so a word that is nothing but ${NAME} goes (env_readings).
    None where env refuses the string and runs nothing: a quote that never
    closes, a $ that opens no ${NAME}, a backslash at the end or before a
    character it takes no escape of, and \\c in double quotes. The env of
    FreeBSD and macOS takes -S too, split by the same rules. Each character
    is looked at once and each word joined once, so the time is in
    proportion to the text (F28's review)."""
    words, word, single, double, i = [], None, False, False, 0

    def add(chars):
        nonlocal word
        if word is None:
            word = []
            words.append(word)
        word.append(chars)

    while i < len(text):
        ch, escaped = text[i], text[i + 1:i + 2]
        if (ch == "'" and not double) or (ch == '"' and not single):
            single, double = (not single, double) if ch == "'" else (single, not double)
            add("")
            i += 1
        elif ch in SPLIT_SPACES and not (single or double):
            word, i = None, i + 1
        elif ch == "#" and word is None:
            break
        elif ch == "\\" and not (single and escaped not in ("\\", "'")):
            if escaped == "c" and not double:
                break
            if escaped == "_" and not double:
                word = None
            elif escaped == "_" or escaped in SPLIT_ESCAPES:
                add(" " if escaped == "_" else SPLIT_ESCAPES[escaped])
            else:
                return None
            i += 2
        elif ch == "$" and not single:
            variable = SPLIT_VARIABLE.match(text, i)
            if variable is None:
                return None
            if not (unset is True or (unset and variable.group()[2:-1] in unset)):
                add(variable.group())
            i = variable.end()
        else:
            add(ch)
            i += 1
    return None if single or double else ["".join(parts) for parts in words]


def split_string_option(word, following):
    """What env's option word gives -S or --split-string, or a start of it,
    to split: the string, attached or following, the next word (None where
    there is none), how many words the option takes, and the letters of a
    bundle before its S, as -i of -iS, which stay as options; None where the
    word gives -S nothing."""
    env = WRAPPERS["env"]
    if word.startswith("--"):
        if long_option_of(word, env) != "split-string":
            return None
        return (word.partition("=")[2], 1, []) if "=" in word else (following, 2, [])
    for at, letter in enumerate(word[1:], start=1):
        if letter in env.short_values:
            if letter != "S":
                return None
            before = [word[:at]] if at > 1 else []
            return (word[at + 1:], 1, before) if word[at + 1:] else (following, 2, before)
    return None


# Which ${NAME}s of env's -S strings a reading takes as unset
# (split_string, readings_of): none, as written, by default.
UNSET = frozenset()


def spliced(rest, read, word, depth, budget):
    """Whether env's option word, at depth, the number of -S strings it came
    out of, gives a string to -S (split_string_option) that is split, its
    words then put on rest, to be read next, as options again, as GNU env
    puts them in the option's place before it reads on (F28: "env -S 'git
    push origin main'" runs the push, and the guard read the string as one
    word and passed it), and the letters of its bundle before the S on
    read. A string env refuses, which runs nothing, is left for the caller
    to read as one word, -S's value, as before the fix. Each word is moved
    once, so many strings cost no more than their words (F28's review: the
    whole command was copied for each). A string -S takes from the next
    word comes out of that word, so its words are one deeper than it, not
    than the -S (the second review: seventeen -S side by side read as
    nested). One nested deeper than MAX_NESTING, or one that takes the
    characters split past budget, a list of what is left, raises TooDeep
    (SPLIT_TOO_DEEP)."""
    given = split_string_option(word, rest[-1][0] if rest else None)
    if given is None or given[0] is None:
        return False
    text, used, (unset, before) = given[0], given[1], (budget[1], given[2])
    words = split_string(text, unset)
    if words is None:
        return False
    if used == 2:
        depth = rest.pop()[1]
    budget[0] -= len(text)
    if depth >= MAX_NESTING or budget[0] < 0:
        raise TooDeep(SPLIT_TOO_DEEP)
    read.extend(before)
    rest.extend((split, depth + 1) for split in reversed(words))
    return True


def env_assigns(word):
    """Whether env reads a word after its options as an assignment to the
    command's environment: any word that holds a =, as GNU env.c and the
    env of FreeBSD and macOS read it, where the shell takes only a name
    and a = (ASSIGNMENT) and runs any other word as a command (F28: "env
    1x=2 git push origin main" runs the push, and the guard read 1x=2 as
    the command)."""
    return "=" in word


def past_lone_dash(rest, read):
    """Move the lone - GNU env reads after its options, a -- among them, as
    -i (env.c) from rest to read, so "env -- - git push" and "env - git
    push" run the push (F28)."""
    if rest and rest[-1][0] == "-":
        read.append(rest.pop()[0])


def read_env_assignments(rest, read, assignments):
    """Move what env reads after its options and before its command from
    rest to read: a lone - (past_lone_dash) and every word env reads as an
    assignment (env_assigns), each added to assignments."""
    past_lone_dash(rest, read)
    while rest and env_assigns(rest[-1][0]):
        in_time()
        assignments.append(rest.pop()[0])
        read.append(assignments[-1])


def describes_only(options):
    """Whether command's options say it runs nothing (DESCRIBES_ONLY), one
    of them alone or in a bundle such as -pv."""
    return any(DESCRIBES_ONLY & set(word[1:]) for word in options if not word.startswith("--"))


class Launch(NamedTuple):
    """What stands before the command itself in a simple command: the
    assignments to its environment, the shell's and env's, each wrapper
    that runs it, by name with its option words, where the command starts,
    None when nothing is run, as after a lone assignment or command -v,
    the words as the wrappers read them, which start indexes: argv itself,
    but where env's -S splits a string into words (F28); each wrapper's
    word as written, /usr/bin/time where its name is time; and why the
    words were not read, where they were not: then nothing is read to run
    (F28's second review: env -S strings nested too deep raised out of the
    reading, and the template's hook, which calls push_arguments itself,
    took the exception for no objection)."""

    assignments: list
    wrappers: list
    start: object
    argv: list
    wrapper_words: tuple = ()
    unread: str = ""
    split: int = 0


def launch(argv, unset=UNSET):
    """argv read as the shell and the wrappers read it, past leading
    assignments, the RESERVED words, function and its name, and the
    WRAPPERS with their options (read_options), env's -S strings split into
    the words env reads (spliced), the ${NAME}s unset names taken as
    unset, and its assignments after them (read_env_assignments) (Launch);
    one read too deep or too wide in its -S strings is read as nothing,
    with why (Launch.unread). The words still to read are kept with the
    next last, so each is moved once, whatever -S puts back; and within a
    run of the hook (refused_push) each reading is kept for the next ask of
    the same words (LAUNCHES; F28's second review: the readings, the moves
    and the places read each command again)."""
    key, kept = (tuple(argv), frozenset(unset)), LAUNCHES.get()
    if kept is not None and key in kept:
        return kept[key]
    try:
        found = read_launch(list(argv), key[1])
    except TooDeep as problem:
        found = Launch([], [], None, list(argv), (), problem.args[0] if problem.args else TOO_DEEP)
    if kept is not None:
        kept[key] = found
    return found


def coproc_named(rest):
    """Whether the word after coproc, the next of rest (read_launch's, the
    next last), is the name it gives the coprocess, which it takes only
    before a compound command, as in coproc X { git push ...; }; before a
    simple command that word is the command, as git is in coproc git push
    (F29's first review)."""
    return len(rest) > 1 and bool(PARAMETER_NAME.fullmatch(rest[-1][0])) and rest[-2][0] in COMPOUND_OPENERS


# The words before a command that may take a name of their own after them,
# which run nothing themselves (read_launch).
NAMING_WORDS = ("function", "coproc")
# The words that open a compound command, after which coproc's name stands.
COMPOUND_OPENERS = {"{", "if", "while", "until", "for", "case", "select", "[[", "(("}


# The readings launch has made in the run under way, by words and the
# names taken as unset, handed out again unchanged; None outside a run.
LAUNCHES = contextvars.ContextVar("LAUNCHES", default=None)


def read_launch(argv, unset):
    """launch's reading, which raises TooDeep where the -S strings nest too
    deep or hold too much to split (spliced)."""
    rest, read = [(word, 0) for word in reversed(argv)], []
    assignments, wrappers, words, budget = [], [], [], [SPLIT_BUDGET, unset]

    def unread():
        return read + [word for word, _ in reversed(rest)]

    def found(start, words_read):
        return Launch(assignments, wrappers, start, words_read, tuple(words), split=SPLIT_BUDGET - budget[0])

    while rest:
        in_time()
        word = rest[-1][0]
        name = program(word)
        if word in RESERVED or ASSIGNMENT.match(word) or word in NAMING_WORDS:
            read.append(rest.pop()[0])
            if ASSIGNMENT.match(word):
                assignments.append(word)
            elif word == "function" and rest or word == "coproc" and coproc_named(rest):
                read.append(rest.pop()[0])
        elif name in WRAPPERS:
            read.append(rest.pop()[0])
            begin = len(read)
            read_options(rest, read, WRAPPERS[name], splits=budget if name == "env" else None)
            for _ in range(WRAPPERS[name].operands):
                if rest:
                    read.append(rest.pop()[0])
            options = read[begin:]
            if name == "command" and describes_only(options):
                return found(None, unread())
            wrappers.append((name, options))
            words.append(word)
            if name == "env":
                read_env_assignments(rest, read, assignments)
        else:
            return found(len(read), unread())
    return found(None, read)


def names_in(argv):
    """The distinct ${NAME}s the words of argv hold, in order."""
    return list(dict.fromkeys(name for word in argv if "${" in word for name in SPLIT_VARIABLE.findall(word)))


def unset_subsets(names):
    """Every set of names a reading takes as unset, the empty one, every
    ${NAME} as written, first: env drops each that is not set, and which are
    is the command's environment's (F28's second review: "env -S '${NOPE}
    -C ${HOME} git push origin main'" runs the push from $HOME where NOPE
    alone is not set)."""
    return [frozenset(chosen) for size in range(len(names) + 1)
            for chosen in itertools.combinations([name[2:-1] for name in names], size)]


# What env's -S may take out of a word between the letters of push, so
# that a word it splits spells push where its text does not (holds_push).
SPLIT_QUOTING = re.compile(r"[\\'\"]")


def holds_push(argv):
    """Whether a word of argv may spell push once env has split it, with
    any of its ${NAME}s unset, at any depth of -S strings: env's -S takes
    out of a word only its quotes, its backslashes and the ${NAME}s that
    are not set, and every escape gives a character that is no letter, or
    a space between words, so a word env may read as push holds push once
    those three are taken out of its text (F28's third review: pu""sh and
    pu${NOPE}sh, inside -S strings, spell push where the text does not)."""
    return any("push" in SPLIT_VARIABLE.sub("", SPLIT_QUOTING.sub("", SPLIT_VARIABLE.sub("", word)))
               for word in argv)


# --- the words bash builds ------------------------------------------------------

# What bash builds a part of a word from once its quotes are gone (F29): a
# parameter, $NAME, ${ }, $1, $@ and the rest of its specials, a
# substitution (SUBSTITUTED), and a pattern's *, ? and [ ], which the name
# of a file fills, and an extended pattern's parenthesis, as in @(main),
# which a word holds only where bash reads it so (pattern_reading); and
# brace expansion, {a,b} and {a..z}, which makes one word many
# (brace_words).
PARAMETER_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
SPECIAL_PARAMETERS = "0123456789@*#?$!-"
BUILT = frozenset("$*?[(" + SUBSTITUTED)


class Hole(NamedTuple):
    """A part of a word bash fills (word_pieces): any text, or one
    character, as a pattern's ?; whether it is an expansion, which may
    give nothing, or a pattern; and, for a bracketed class, what stands
    between its brackets, which says the characters it may give
    (hole_allows), None for any."""

    any: bool
    expansion: bool
    allows: str = None


EXPANSION_HOLE, PATTERN_HOLE, CHARACTER_HOLE = Hole(True, True), Hole(True, False), Hole(False, False)


def expansion_end(word, i):
    """Where the expansion bash makes of the $ at word[i] ends, a ${ } at
    the } that closes it, counting the braces inside, or at the end of the
    word where none does; None where the $ is a character of its own, as it
    is before anything but a name, a { or a special parameter."""
    following = word[i + 1:i + 2]
    if following == "{":
        depth = 0
        for brace in BRACES.finditer(word, i + 1):
            depth += 1 if brace.group(0) == "{" else -1
            if not depth:
                return brace.end()
        return len(word)
    name = PARAMETER_NAME.match(word, i + 1)
    if name:
        return name.end()
    return i + 2 if following and following in SPECIAL_PARAMETERS else None


BRACES = re.compile(r"[{}]")


def word_pieces(word):
    """A word as bash builds it, its quotes gone: the runs of its own
    characters, as text, and in their places the holes bash fills (Hole):
    each expansion, any text, and each part of a pattern, * any text and
    ? or a bracketed class one character. The guard keeps no quotes, so a
    $ or a * that was quoted reads as one bash fills, which refuses too
    much at worst. The last few words' pieces are kept (pieces_of), since
    each rule that judges a word asks for them again, and a word of
    MAX_COMMAND characters took a second each time (F29's second
    review)."""
    return list(pieces_of(word))


@functools.lru_cache(maxsize=16)
def pieces_of(word):
    """word_pieces's reading of a word, as a tuple."""
    pieces, text, i, no_class_after, closes = [], [], 0, len(word), None

    def hole(kind):
        if text:
            pieces.append("".join(text))
            text.clear()
        pieces.append(kind)

    while i < len(word):
        ch = word[i]
        end = expansion_end(word, i) if ch == "$" else None
        close = word.find("]", i + 2) if ch == "[" and i < no_class_after else -1
        if ch in "@?*+!" and closes is None and "(" in word:
            closes = closing_parentheses(word, escapes=False)
        group = pattern_end(word, i, closes) if ch in "@?*+!" and closes else None
        if end is not None or ch == SUBSTITUTED:
            hole(EXPANSION_HOLE)
            i = end or i + 1
        elif group:
            hole(PATTERN_HOLE)
            i = group
        elif ch in "*?":
            hole(PATTERN_HOLE if ch == "*" else CHARACTER_HOLE)
            i += 1
        elif close != -1:
            hole(Hole(False, False, word[i + 1:close]))
            i = close + 1
        else:
            no_class_after = i if ch == "[" else no_class_after
            text.append(ch)
            i += 1
    if text:
        pieces.append("".join(text))
    return tuple(pieces)


def pattern_end(word, i, closes):
    """Where the extended pattern that opens at word[i], as @(main|x), ends,
    after the parenthesis that closes it (closes, the word's
    closing_parentheses); None where none opens there or none closes
    (word_pieces)."""
    if not word.startswith("(", i + 1) or word[i - 1:i] == "$":
        return None
    return closes.get(i + 1)


def hole_allows(hole, ch):
    """Whether a hole may give the character ch: any but a bracketed
    class, which gives what its brackets hold, as [+] gives + and [!a] any
    but a, a file's name filling it (F29's second review). A class that
    names a character class, as [[:punct:]], is read as allowing any,
    which refuses too much at worst."""
    if hole.allows is None or "[:" in hole.allows:
        return True
    body = "!" + hole.allows[1:] if hole.allows.startswith("^") else hole.allows
    return fnmatch.fnmatchcase(ch, "[" + body + "]")


def is_built(word):
    """Whether bash builds a part of the word: it holds a character that
    may open a part bash builds (BUILT), and one does (word_pieces), so a
    $ that ends a word, as the whole reading leaves of --no-verify$(true),
    is the character it is."""
    return not BUILT.isdisjoint(word) and any(isinstance(piece, Hole) for piece in word_pieces(word))


def may_be(word, text):
    """Whether bash may build the word into text: its own characters, in
    order, around holes each filled as it may be (word_pieces), as
    git pu$(true)sh and git $(echo push) run git push, and a push to
    ma$(true)in, m${V:-ai}n or mai? one to main (F29). Each hole may give
    any text but a bracketed class, which gives one of the characters it
    holds (hole_allows); a hole is never read as a space between words,
    which a variable's value or a command's output may hold too. Where a
    hole of a pattern may give a + or a colon, as a file's name may hold
    one, the refspec it makes is judged by may_open_with_plus and
    colon_given."""
    if not is_built(word):
        return word == text
    return pieces_may_be(word_pieces(word), text)


def pieces_may_be(pieces, text):
    """may_be's reading of a word's pieces (word_pieces)."""
    if not pieces:
        return not text
    if sum(len(piece) for piece in pieces if isinstance(piece, str)) > len(text):
        return False
    places = {0}
    for piece in pieces:
        if isinstance(piece, str):
            places = {at + len(piece) for at in places if text.startswith(piece, at)}
        elif piece.any:
            places = set(range(min(places), len(text) + 1))
        else:
            places = {at + 1 for at in places if at < len(text) and hole_allows(piece, text[at])}
        if not places:
            return False
    return len(text) in places


def expansions_alone(word):
    """Whether the word is made of expansions alone, with no character of
    its own and no pattern, as "$(git branch --show-current)", "$BRANCH"
    and "${PREFIX}${N}" are: a value the guard cannot know, read as a
    variable's is, and, where bash may make it nothing, as nothing
    (unset_value). Two expansions side by side may make main as one may
    (F29's first review), so the guard reads them as it reads one, and
    refuses only a name that holds a character of its own beside them, as
    ma$(true)in does, or a pattern, which the names of files fill."""
    return is_built(word) and all(piece == EXPANSION_HOLE for piece in word_pieces(word))


# A ${ } as unset_value reads its start: an indirection ! or a length #,
# the parameter, an array's subscript after it, and the operator after
# that, of which -, :-, = and := give their word where the parameter is
# unset, as ${NOPE[@]:-main} gives main (F29's first review), and + and :+
# where it is set (defaults_given). The subscript stops at a }, a newline
# or its bound, MOST_SUBSCRIPT characters: read to the next ] however far,
# ${A[ written again and again was read in a time that grew as the square
# of the command, 10.8 seconds at MAX_COMMAND (F29's second review).
MOST_SUBSCRIPT = 64
BRACED_HEAD = re.compile(r"\$\{([!#]?)([A-Za-z_][A-Za-z0-9_]*(?:\[[^\]}\n]{0,%d}\])?|[0-9]+|[@*#?$!-])?(:?[-=+?])?"
                         % MOST_SUBSCRIPT)
GIVES_ITS_WORD = ("-", ":-", "=", ":=")
GIVES_ITS_WORD_SET = ("+", ":+")


def unset_value(word, carried=None):
    """The word as bash builds it where every variable is unset and every
    substitution prints nothing (F29): each $NAME, special parameter and
    substitution gone, and each ${ } replaced by what it gives then, the
    word of a ${NAME:-word} and its like, read the same way, a $'...' in it
    as the text it decodes to (ansi_c_text), which bash gives there even
    inside double quotes, and nothing for any other, so ${NOPE:-push} is
    push and pu${NOPE}sh push. A world bash may run in, so a push found in
    it is one bash may run. carried, where given, is the ${ }s a word
    before this one left open, which this one goes on with and leaves as
    it finds them (unset_words), each with whether it gives its word, the
    braces inside it and how many of those open around it give none."""
    out, stack, i = [], [] if carried is None else carried, 0
    dropping = stack[-1][2] if stack else 0
    while i < len(word):
        ch = word[i]
        if ch == "$" and word.startswith("${", i):
            head = BRACED_HEAD.match(word, i)
            keeps = not head.group(1) and head.group(3) in GIVES_ITS_WORD
            dropping += not keeps
            stack.append([keeps, 0, dropping])
            i = head.end()
            continue
        if stack and (decoded := ansi_c_in_braces(word, i)):
            out += [] if dropping else [decoded[0]]
            i = decoded[1]
            continue
        end = expansion_end(word, i) if ch == "$" else None
        if end is not None or ch == SUBSTITUTED:
            i = end or i + 1
            continue
        if stack and ch == "}" and not stack[-1][1]:
            stack.pop()
            dropping = stack[-1][2] if stack else 0
            i += 1
            continue
        if stack and ch in "{}":
            stack[-1][1] += 1 if ch == "{" else -1
        if not dropping:
            out.append(ch)
        i += 1
    return "".join(out)


def defaults_given(word, carried=None):
    """The word as bash builds it where each ${NAME<op>word} whose operator
    gives its word, -, :-, = and := where the parameter is unset and + and
    :+ where it is set, gives it, a $'...' in it as the text it decodes to,
    and every other expansion stands as written, a hole (F29's second
    review). A default is the command's own text, so ${A:-ma}$(echo in),
    which bash makes main where A is unset, reads as ma$(echo in) does,
    and ${A:-+}$B as +$B. carried, where given, is the ${ }s a word before
    this one left open, each the count of braces open inside it, which
    this one goes on with (defaults_way)."""
    out, stack, i = [], [] if carried is None else carried, 0
    while i < len(word):
        ch = word[i]
        if word.startswith("${", i):
            head = BRACED_HEAD.match(word, i)
            if not head.group(1) and head.group(3) in GIVES_ITS_WORD + GIVES_ITS_WORD_SET:
                stack.append(0)
                i = head.end()
                continue
            end = expansion_end(word, i)
            out.append(word[i:end])
            i = end
            continue
        if stack and (decoded := ansi_c_in_braces(word, i)):
            out.append(decoded[0])
            i = decoded[1]
            continue
        if stack and ch == "}" and not stack[-1]:
            stack.pop()
            i += 1
            continue
        if stack and ch in "{}":
            stack[-1] += 1 if ch == "{" else -1
        out.append(ch)
        i += 1
    return "".join(out)


def defaults_way(way):
    """A simple command's words with each default given (defaults_given),
    a ${ } left open at a word's end going on into the words after it, as
    unset_words reads it, and a word that is then nothing dropped; the
    words themselves where none holds a ${ }."""
    if not any("${" in word for word in way):
        return list(way)
    carried, words = [], []
    for word in way:
        in_time()
        value = defaults_given(word, carried)
        if value or not word:
            words.append(value)
    return words


def kept_as_written(word, value):
    """Whether a word bash builds, whose unset value is only colons and
    pluses, is read as written in the unset world (unset_words): only one
    that opens and ends with a hole, its defaults given (defaults_given),
    as $S:$D does, so :$X, $D: and ${A:-:}$B are read as the : bash makes
    of them where the hole is nothing (F29's second review)."""
    if not value or value.strip(":+") or not is_built(word):
        return False
    pieces = word_pieces(defaults_given(word))
    return isinstance(pieces[0], Hole) and isinstance(pieces[-1], Hole)


def ansi_c_in_braces(word, i):
    """The $'...' that opens at word[i], inside a ${ }, as the text it
    decodes to (ansi_c_text) and where it ends; None where none opens there
    or none closes. bash decodes one there even inside double quotes, where
    plain_text keeps it as written (unset_value)."""
    if not word.startswith("$'", i):
        return None
    close = ansi_c_close(word, i + 2)
    return (ansi_c_text(word[i + 2:close - 1]), close) if close else None


def unset_words(way):
    """A simple command's words where every variable is unset and every
    substitution prints nothing (unset_value), a word that is then nothing
    dropped. A ${ } that a blank stands in, unquoted, is split there into
    words by shlex as bash splits what it gives, so the words after it go
    on with it until it closes: ${NOPE:-git push origin main} runs git push
    origin main, while its quoted form stays one word, which runs no git
    (F29's first review). And a word with holes that is then only colons and
    pluses, as $B:$B, is kept as written, as one whole expansion is: read
    as nothing in each hole it is the refspec :, which pushes every branch
    both sides have, but the session writes it for one branch whose name it
    reads, as git push origin "$(git branch --show-current)" does (F29's
    first review), and a + in it, as in +$B, still forces."""
    carried, words = [], []
    for word in way:
        in_time()
        value = unset_value(word, carried)
        if kept_as_written(word, value):
            value = word
        if value or not word:
            words.append(value)
    return words


# How many words brace expansion may add to one simple command's before the
# guard stops expanding it: past it a command that may run git is refused
# unread (TOO_MANY_WORDS), since a push may stand in a word not made yet.
MOST_WORDS = 256
TOO_MANY_WORDS = (f"a command whose braces expand to more than {MOST_WORDS} words, or to more text than the guard "
                  "reads, is refused when it may run git: the guard cannot read every word in time; write the words "
                  "out (.claude/rules/board.md)")
SEQUENCE = re.compile(r"(-?[0-9]+)\.\.(-?[0-9]+)(?:\.\.(-?[0-9]+))?|([A-Za-z])\.\.([A-Za-z])(?:\.\.(-?[0-9]+))?")


class TooMany(Exception):
    """Brace expansion that makes more words or text than the guard reads
    (MOST_WORDS, SPLIT_BUDGET), or nests deeper than MAX_NESTING."""


def brace_pairs(word):
    """The braces of a word that brace expansion may expand, as (start,
    end, commas), in the order they open: each { with the } that closes it,
    the commas directly inside it, and none that is a ${ } or stands in one,
    which bash expands as a parameter (x${a,b}{c,d} makes xc and xd). Read
    in one pass, each brace matched by a stack, as bash's brace_gobbler
    matches them from each {."""
    stack, pairs, parameters = [], [], 0
    for mark in BRACE_MARKS.finditer(word):
        i, ch = mark.start(), mark.group(0)
        if ch == "{":
            parameter = i > 0 and word[i - 1] == "$"
            stack.append((i, [], parameter))
            parameters += parameter
        elif ch == "}" and stack:
            start, commas, parameter = stack.pop()
            parameters -= parameter
            if not parameter and not parameters:
                pairs.append((start, i + 1, commas))
        elif ch == "," and stack:
            stack[-1][1].append(i)
    return sorted(pairs)


BRACE_MARKS = re.compile(r"[{},]")


def alternatives_of(word, start, end, commas):
    """What a brace pair gives in its place: the span between each two of
    its commas, as (start, end), each expanded in turn, or a sequence's
    terms as text (sequence_terms); None where it gives nothing, as {a}
    and {1..a} do, which bash keeps as written."""
    if commas:
        cuts = [start, *commas, end - 1]
        return [(cuts[k] + 1, cuts[k + 1]) for k in range(len(cuts) - 1)]
    sequence = SEQUENCE.fullmatch(word, start + 1, end - 1)
    return sequence_terms(sequence) if sequence else None


def sequence_terms(sequence):
    """The terms of a brace sequence, {1..10..3} or {a..e}, as bash makes
    them: from the first to the last by the step, its sign ignored, numbers
    padded with zeros to the longer of the two where either starts with a
    zero; TooMany where they are more than MOST_WORDS."""
    first, last, step, letter_first, letter_last, letter_step = sequence.groups()
    low, high = (int(first), int(last)) if first is not None else (ord(letter_first), ord(letter_last))
    step = abs(int((step if first is not None else letter_step) or 1)) or 1
    if abs(high - low) // step + 1 > MOST_WORDS:
        raise TooMany()
    terms = range(low, high + (1 if high >= low else -1), step if high >= low else -step)
    if first is None:
        return [chr(term) for term in terms]
    padded = any(len(number.lstrip("-")) > 1 and number.lstrip("-").startswith("0") for number in (first, last))
    width = max(len(first), len(last)) if padded else 0
    return [f"{term:0{width}d}" for term in terms]


def brace_words(word):
    """The words bash's brace expansion makes of a word, left to right, as
    its brace_expand makes them: each leftmost brace pair it expands
    (brace_pairs, alternatives_of) giving each of its alternatives,
    themselves expanded, between the text before it and the words of the
    text after it, so p{u,}sh makes push and psh and git {push,-v} runs git
    push -v (F29). Each brace is matched once and each pair read once,
    however deep they nest. TooMany past MOST_WORDS words, past SPLIT_BUDGET
    characters, or nested past MAX_NESTING."""
    if "{" not in word:
        return [word]
    pairs = brace_pairs(word)
    return expanded_span(word, pairs, [pair[0] for pair in pairs], 0, len(word), 0)


def expanded_span(word, pairs, starts, low, high, depth):
    """The words brace expansion makes of word[low:high] (brace_words),
    pairs and their starts the word's brace pairs, at depth."""
    if depth > MAX_NESTING:
        raise TooMany()
    made, at, k = [""], low, bisect.bisect_left(starts, low)
    while k < len(pairs) and pairs[k][0] < high:
        in_time()
        start, end, commas = pairs[k]
        alternatives = alternatives_of(word, start, end, commas) if end <= high else None
        if alternatives is None:
            k += 1
            continue
        given = [made_word for alternative in alternatives for made_word in (
            [alternative] if isinstance(alternative, str) else
            expanded_span(word, pairs, starts, alternative[0], alternative[1], depth + 1))]
        made = [before + word[at:start] + after for before in made for after in given]
        if len(made) > MOST_WORDS or sum(map(len, made)) > SPLIT_BUDGET:
            raise TooMany()
        at, k = end, bisect.bisect_left(starts, end)
    return [before + word[at:high] for before in made]


def word_ways(argv):
    """The ways bash may build the words of a simple command, its quotes
    gone (F29): as written, each word whole, the reading every other part of
    the guard makes; with its braces expanded (brace_words), where it holds
    any; and each of those where every variable is unset and every
    substitution prints nothing (unset_value), a word that is then nothing
    dropped, as bash drops an unquoted word that expands to nothing, so a
    command whose first word is $NOPE or $(true) runs the git after it; and
    whether the braces made more than the guard reads (TooMany). Within a
    run of the hook each argv's ways are kept for the next ask
    (BUILT_WAYS), as launch keeps its readings."""
    key, kept = tuple(argv), BUILT_WAYS.get()
    if kept is not None and key in kept:
        return kept[key]
    found = built_ways(argv)
    if kept is not None:
        kept[key] = found
    return found


# The ways word_ways has found in the run under way, by words, and the
# readings readings_of has made, handed out again unchanged; None outside a
# run (refused_push).
BUILT_WAYS = contextvars.ContextVar("BUILT_WAYS", default=None)
READINGS = contextvars.ContextVar("READINGS", default=None)


def built_ways(argv):
    """word_ways's reading of argv, unkept."""
    ways, too_many, text = [list(argv)], False, "".join(argv)
    if BUILT.isdisjoint(text) and "{" not in text:
        return ways, too_many
    if "{" in text:
        try:
            braced = [made for word in argv for made in brace_words(word) if made or not word]
            if len(braced) > len(argv) + MOST_WORDS:
                raise TooMany()
        except TooMany:
            too_many = True
        else:
            ways += [braced] if braced != ways[0] else []
    for way in list(ways):
        given = defaults_way(way)
        ways += [given] if given not in ways else []
    for way in list(ways):
        if any(not BUILT.isdisjoint(word) for word in way):
            unset = unset_words(way)
            ways += [unset] if unset not in ways else []
    return ways, too_many


def names_git(word):
    """Whether a command's first word may run git: git by any directory, or
    a word bash may build into git (may_be), as $(which git) or gi?."""
    return program(word) == "git" or (is_built(word) and may_be(program(word), "git"))


def readings_of(argv):
    """The ways a simple command may run (Launch): in each way bash may
    build its words (word_ways, F29), each of the ways env may run it then
    (env_ways), each reading once; and, where its braces make more words
    than the guard reads and its first word may run git (names_git), a way
    read as nothing, with why (TOO_MANY_WORDS). Within a run of the hook
    each argv's readings are kept for the next ask (READINGS)."""
    key, kept = tuple(argv), READINGS.get()
    if kept is not None and key in kept:
        return kept[key]
    found = built_readings(argv)
    if kept is not None:
        kept[key] = found
    return found


def built_readings(argv):
    """readings_of's reading of argv, unkept."""
    ways, too_many = word_ways(argv)
    if len(ways) == 1 and not too_many:
        return env_ways(ways[0])
    found = {}
    for way in ways:
        for reading in env_ways(way):
            found.setdefault((tuple(reading.argv), reading.start, reading.unread), reading)
    readings = list(found.values())
    if too_many and may_push_braced(readings[0]):
        readings.append(readings[0]._replace(start=None, unread=TOO_MANY_WORDS))
    return readings


def may_push_braced(found):
    """Whether a reading (Launch) whose braces were not expanded may run a
    push once they are: its first word may run git (names_git) and its
    verb, every brace pair in it a hole (braced_as_holes), may be push, as
    the verb {push,-v} may be, where {commit,-v} may not."""
    verb = verb_at(found)
    return verb is not None and may_be(braced_as_holes(found.argv[verb]), "push")


def braced_as_holes(word):
    """The word with each brace pair that brace expansion may expand
    (brace_pairs), outermost, standing as one hole (SUBSTITUTED), since what
    it gives is not read."""
    if "{" not in word:
        return word
    kept, at = [], 0
    for start, end, _ in brace_pairs(word):
        if start >= at:
            kept += [word[at:start], SUBSTITUTED]
            at = end
    return "".join(kept) + word[at:]


def env_ways(argv):
    """The ways env may run a simple command (Launch): as launch reads it,
    each ${NAME} in an -S string kept as written; and, where the strings
    hold any, once more for each set of them unset (unset_subsets), env
    expanding an unset one to nothing, so a word that is nothing but
    ${NAME} goes and one that starts with it starts with what follows
    (F28's review: "env -S '${NOPE} git push origin main'" runs the push
    where NOPE is not set). Every mix is read, since which words env drops
    decides which option takes which word as its value, so a push only one
    mix runs, as "env -S '${A}-u ${B} ${C}-u ${D}' x git push origin main"
    does where B alone is set, is found by no other: reading the names all
    kept and all unset alone would let it through (F28's third review).
    Past MOST_NAMES distinct names, or where the mixes would split more
    than SPLIT_BUDGET characters in all, the mixes are not read, and a
    command whose words may spell push (holds_push) is read as nothing,
    with why (SPLIT_TOO_MANY_NAMES), since a push to judge may stand in a
    mix the guard has not read; one whose words cannot spell push runs no
    push in any mix, and is read as written."""
    kept = launch(argv)
    names = names_in(argv)
    if not names or kept.argv == list(argv) or kept.unread:
        return [kept]
    if len(names) > MOST_NAMES or 2 ** len(names) * max(kept.split, 1) > SPLIT_BUDGET:
        return [kept] + ([kept._replace(start=None, unread=SPLIT_TOO_MANY_NAMES)] if holds_push(argv) else [])
    found, spent = {}, 0
    for unset in unset_subsets(names):
        reading = launch(argv, unset)
        spent += reading.split
        if spent > SPLIT_BUDGET:
            return [kept, kept._replace(start=None, unread=SPLIT_TOO_DEEP)]
        found.setdefault(tuple(reading.argv), reading)
    return list(found.values())


def env_readings(argv):
    """The words of each way env may run a simple command (readings_of),
    argv itself for a way that could not be read, whose push push_arguments
    then refuses unread; the guard judges every push any of them finds."""
    return [list(argv) if reading.unread else reading.argv for reading in readings_of(argv)]


def command_start(argv):
    """Where the command itself starts in the words as the wrappers read
    them, argv but where env's -S splits a string (Launch.argv), or None
    when nothing is run (launch)."""
    return launch(argv).start


class Unreadable(str):
    """A word push_arguments gives back in place of a push's arguments
    where it could not read the command, which holds a push: its text is
    why, which refusal gives back, so the template's hook, which calls
    push_arguments and refusal itself, refuses the command too (F28's
    second review)."""


class Readings(list):
    """A push's arguments as one way of reading the command finds them, the
    first, with others, those every other way that finds a push finds,
    each of which refusal judges (F28's second review: the template's hook
    judged only one reading)."""

    def __init__(self, first, others):
        super().__init__(first)
        self.others = others


def push_arguments(argv):
    """The arguments after "push" when argv is a git push, else None, in
    every way env may run it (readings_of): the one way's alone, or, where
    ways differ, Readings, all of them; and [Unreadable(why)] where a way
    could not be read, whatever its words, since the words env would run
    are what was not read (F28's third review: a push spelled pu""sh in an
    -S string nested too deep passed, its text holding no push). The hook
    asks only of a command that holds push somewhere in its text, so one
    that cannot push is not refused for this. git may stand behind
    assignments and wrappers (command_start), and a query of git's own, as
    --exec-path alone, runs no command at all."""
    readings = readings_of(argv)
    unread = next((reading.unread for reading in readings if reading.unread), "")
    if unread:
        return [Unreadable(unread)]
    found = list(dict.fromkeys(tuple(args) for args in map(arguments_of, readings) if args is not None))
    if not found:
        return None
    return list(found[0]) if len(found) == 1 else Readings(found[0], [list(args) for args in found[1:]])


def arguments_of(found):
    """The arguments after "push" in one reading (Launch), else None. A
    command's first word bash may build into git, and a verb it may build
    into push, are read as git and push (names_git, may_be; F29: git
    $(echo push) origin main and $(which git) push origin main run the
    push), so a verb such as pu$(x) is read as push, the words after it as
    its arguments."""
    verb = verb_at(found)
    word = None if verb is None else found.argv[verb]
    return found.argv[verb + 1:] if word is not None and (word == "push" or may_be(word, "push")) else None


def verb_at(found):
    """Where git's verb stands in one reading (Launch), past git's own
    options, each with its value where it takes one (GIT_VALUE_OPTIONS),
    where its first word may run git (names_git); None where none does, or
    a query of git's own runs no command (GIT_QUERIES). After a first word
    that is not git as written but may be built into it, an option git does
    not know is one git refuses, running nothing (GIT_FLAGS): read with no
    such care, env -S '${NOPE} -u ${HOME} git push ...', where env reads
    each ${NAME} as set, ran git -u with the value of HOME as its verb, a
    world in which git runs nothing."""
    start, argv = found.start, found.argv
    if start is None or not names_git(argv[start]):
        return None
    as_written, i = program(argv[start]) == "git", start + 1
    while i < len(argv):
        word = argv[i]
        if word in GIT_QUERIES or word.startswith("--list-cmds"):
            return None
        if word in GIT_VALUE_OPTIONS:
            i += 2
        elif word.startswith("-"):
            if not (as_written or word in GIT_FLAGS or word.startswith(GIT_ATTACHED)):
                return None
            i += 1
        else:
            return i
    return None


def assigns_repository(word):
    """Whether an assignment, or a name export hands on, is to a variable
    that names the repository (GIT_REPOSITORY_VARIABLES)."""
    return word.split("=", 1)[0].rstrip("+") in GIT_REPOSITORY_VARIABLES


def joined(place, path, physical=False):
    """The directory a path names from place: an absolute path in place of
    it, whatever place is, a relative one joined on, ~ the home directory,
    and an empty path, which git -C and bash's cd both take as no change,
    place itself; None when the path holds what the shell expands first
    (EXPANDED) or starts with a tilde of bash's that names one of the
    shell's own directories (SHELL_TILDE), when it is relative and place is
    None, and when the system cannot read it, as a ~ before a NUL byte or a
    character no file name can hold, which expanduser refuses with
    ValueError (F21's third review: until then such a path raised, and the
    hook with it, and an absolute path after a place this module could not
    tell was one it could not tell either), and when the path it names is
    too long for the system to take in one call (PATH_MAX). A logical
    path, the one cd keeps, is normalised as the shell normalises it, a ..
    taking back the step before it; a physical one, as git -C and env -C change directory,
    is left as written for the system to resolve, since there a .. after a
    symbolic link climbs from where the link leads."""
    if EXPANDED & set(path) or SHELL_TILDE.match(path):
        return None
    if not path:
        return place
    try:
        path = os.path.expanduser(path)
    except ValueError:
        return None
    if not os.path.isabs(path):
        if place is None:
            return None
        path = os.path.join(place, path)
    path = path if physical else os.path.normpath(path)
    return None if len(path) >= PATH_MAX or len(path.encode("utf-8", "surrogatepass")) >= PATH_MAX else path


def names_chdir(word):
    """Whether one of env's long option words names --chdir, in full or by
    any start of it that names no other option of env's, as --ch does
    (long_option_of)."""
    return long_option_of(word, WRAPPERS["env"]) == "chdir"


def env_directory(options):
    """The directory env's option words name for the command it runs, as
    GNU env reads -C DIR, -CDIR, a bundle that ends in C, and --chdir or a
    start of it (names_chdir, F26: "env --ch=../other" was judged where env
    was run), with its value attached or as the next word: the last one
    given, since env changes directory once, to the last directory it was
    given, from where it was run (F21's third review: "env -C sub -C
    ../other" runs in ../other); None where none is given."""
    found, env, i = None, WRAPPERS["env"], 0
    while i < len(options):
        word = options[i]
        i += 1
        if word == "--":
            break
        if word.startswith("--"):
            value = word.partition("=")[2]
            if long_value_is_next(word, env) and i < len(options):
                value, i = options[i], i + 1
            if names_chdir(word):
                found = value
            continue
        for at, letter in enumerate(word[1:], start=1):
            if letter in env.short_values:
                value = word[at + 1:]
                if not value and i < len(options):
                    value, i = options[i], i + 1
                if letter == "C":
                    found = value
                break
    return found


def chdir_paths(argv):
    """The paths the git push argv changes directory by before it pushes,
    in the order chdir reads them: the directory each env before it names
    (env_directory), then each -C git is given; None where the push names
    its repository by a path, with --git-dir, --bare or a variable set on
    it (GIT_REPOSITORY_VARIABLES), so it goes where this module does not
    look."""
    found = launch(argv)
    argv = found.argv
    if found.start is None or any(assigns_repository(word) for word in found.assignments):
        return None
    named = (env_directory(options) for name, options in found.wrappers if name == "env")
    paths = [directory for directory in named if directory is not None]
    i = found.start + 1
    while i < len(argv) and argv[i] != "push":
        word = argv[i]
        if word.split("=", 1)[0] in GIT_REPOSITORY_OPTIONS:
            return None
        if word == "-C":
            paths.append(argv[i + 1] if i + 1 < len(argv) else "")
        i += 2 if word in GIT_VALUE_OPTIONS else 1
    return paths


def push_place(argv, here):
    """The directory the git push argv runs in, from here, the directory
    its shell runs it in (placed_commands): the paths it changes directory
    by (chdir_paths), joined on in turn as chdir reads them (joined,
    physical); None where joined cannot tell, or where the push names its
    repository by a path."""
    paths = chdir_paths(argv)
    if paths is None:
        return None
    place = here
    for path in paths:
        place = joined(place, path, physical=True)
    return place


# --- where the shell runs each command ------------------------------------------

class Stacked:
    """The directory stack after a pushd: the directory it stacked (top) on
    the stack before it (below), () being the empty stack. A stack is
    made by stacked alone, which makes one object of every stack of the
    same entries, so a pushd, a popd and comparing two stacks each take the
    same time however deep the stack is; a stack kept as a tuple took a
    time in its depth for each, and a run of pushd one in the square of its
    length (F25's third review: 29,000 of "pushd a;" took 68 seconds)."""

    __slots__ = ("top", "below", "__weakref__")

    def __init__(self, top, below):
        self.top, self.below = top, below


# Every stack made and still held, by its top and the stack below it.
STACKS = weakref.WeakValueDictionary()


def stacked(below, top):
    """The stack of top on below (Stacked), one object for each."""
    key = (top, below)
    found = STACKS.get(key)
    if found is None:
        found = STACKS[key] = Stacked(top, below)
    return found


class Shell(NamedTuple):
    """Where a shell stands as it runs a command: its directory and the
    directories pushd stacked (Stacked), each None where this module cannot
    tell; whether a cd, pushd or popd ran before, made or failed, after which
    this module cannot tell that the shell stands where it reads it
    (after_a_move); and whether a variable exported names the
    repository (GIT_REPOSITORY_VARIABLES), which holds wherever the shell
    moves after, since a move changes the directory and not the variable
    (followed_and_unmoved); and whether such a variable may have been set,
    alone, which git sees only where it was exported before, or behind a
    time the guard cannot tell from bash's keyword
    (maybe_naming_repository)."""

    directory: object
    stack: object = ()
    moved: bool = False
    named_repository: bool = False
    maybe_named: bool = False


def operators(run):
    """The shell's operators in a run of them, ;; and && and the like as
    one each: ")", "|" for ")|"."""
    found, i = [], 0
    while i < len(run):
        size = 2 if run[i:i + 2] in ("&&", "||", "|&", ";;", ";&") else 1
        found.append(run[i:i + size])
        i += size
    return found


def operands_of(words):
    """A builtin's operands, past its options and a -- that ends them; a
    lone - is an operand."""
    i = 0
    while i < len(words) and words[i].startswith("-") and words[i] != "-":
        i += 1
        if words[i - 1] == "--":
            break
    return words[i:]


def rotates(word):
    """Whether pushd's or popd's operand is +N or -N, which turns the stack."""
    return word[:1] in "+-" and word[1:].isascii() and word[1:].isdigit()


def reached(shell, target, stack, is_directory):
    """The shell after a move to target: there, with stack, when target is
    a directory as the hook runs (is_directory), or one this module cannot
    tell; where it stood when target is not, since cd and pushd fail there
    and the shell runs what follows a ; or a newline from where it stood,
    and before F21's second review such a push was judged in a directory
    that holds nothing, so a lease from the clone the shell never left
    went unjudged."""
    if target is not None and not is_directory(target):
        return shell
    return Shell(target, stack)


def cd(shell, operands, is_directory):
    """Where cd takes the shell: to its operand, the home directory with
    none, along the logical path the shell keeps (joined), and where
    reached says for a target that is not a directory; where this module
    cannot tell for -, the directory before, which it does not track, and
    with -P, which takes the shell along the physical path, so that a
    later .. climbs from where a symbolic link led. -L, -e and -@ change
    nothing here."""
    words = operands_of(operands)
    options = operands[:len(operands) - len(words)]
    target = words[0] if words else "~"
    if target == "-" or any("P" in word[1:] for word in options):
        return Shell(None, shell.stack)
    return reached(shell, joined(shell.directory, target), shell.stack, is_directory)


def pushd(shell, operands, is_directory):
    """pushd DIR moves the shell to DIR and stacks where it was, where
    reached says for a DIR that is not a directory; with -n it stacks DIR
    and moves nothing, and with no operand, -, +N or -N it turns the
    stack, which this module does not follow."""
    words, moves = list(operands), True
    while words and words[0] == "-n":
        words, moves = words[1:], False
    if words and words[0] == "--":
        words = words[1:]
    if not words or words[0] == "-" or rotates(words[0]):
        return Shell(None, None)
    if not moves:
        return Shell(shell.directory, None)
    stack = None if shell.stack is None else stacked(shell.stack, shell.directory)
    return reached(shell, joined(shell.directory, words[0]), stack, is_directory)


def popd(shell, operands):
    """popd moves the shell back to where pushd last stacked it; with -n,
    +N or -N it takes another entry off, or one without the move, and from
    a stack this command did not fill it moves where this module cannot
    tell."""
    if operands == ["-n"]:
        return Shell(shell.directory, None)
    if operands or not shell.stack:
        return Shell(None, None)
    return Shell(shell.stack.top, shell.stack.below)


# The words before a command that leave it run by the shell itself, so a
# cd behind them moves the shell: command, and time, which bash reads as a
# keyword that times what follows in the same shell.
IN_THE_SHELL = ("command", "time")


def after_move(before, after):
    """The shell after a cd, pushd or popd, made or failed, from the shell
    before it: one that has moved, and whose repository a variable names
    where it named it before, since a move changes the directory and not
    the variable."""
    return after._replace(moved=True, named_repository=before.named_repository, maybe_named=before.maybe_named)


def naming_repository(shell):
    """The shell once a variable that names the repository is set alone or
    exported (GIT_REPOSITORY_VARIABLES): git goes to a repository this
    module does not follow, which puts the directory out of its reading,
    and the variable names it after any later move or -C (Shell)."""
    return shell._replace(directory=None, named_repository=True)


def sets_the_shells_own(found):
    """Whether the assignments of a simple command that runs nothing
    (Launch) set the shell's own variables, as a lone GIT_DIR=x does, and
    after bash's time, which times the assignment in the shell: they do
    not behind any other wrapper, since env GIT_DIR=x, with no command,
    prints its environment and leaves the shell's as it was (F28: until
    the fix such a line was read as GIT_DIR set in the shell, and every
    push after it judged by its words and in the session's own clone
    alone, so a lease by git -C onto another worker's claim that a clone
    elsewhere shows passed). bash's keyword is the word time itself, so
    /usr/bin/time, the program, which runs the assignment as a command,
    sets nothing (F28's review: the guard matched the program's name)."""
    return all(word == "time" for word in found.wrapper_words)


def behind_time_alone(found):
    """Whether a simple command's wrappers are bash's time and nothing else,
    each written time (sets_the_shells_own): the lexer drops quotes, so a
    "time" or \\time, which bash runs as the program and not as its
    keyword, reads the same, and the guard cannot tell whether the
    assignment set the shell's variable (maybe_naming_repository)."""
    return bool(found.wrapper_words) and all(word == "time" for word in found.wrapper_words)


def set_alone(found):
    """Whether a simple command that runs nothing is its assignments alone,
    behind no wrapper: bash sets each as a variable of the shell's, which
    git sees only where the name was exported before, so the guard cannot
    tell whether GIT_DIR=x names the repository of a push after it (F29:
    until the fix a lone GIT_DIR=x was read as naming it for every push
    after, each judged by its words and in the session's own clone alone,
    while git ran a push by -C in the clone the -C names, as it does after
    time GIT_DIR=x, which the guard judged there too)."""
    return not found.wrapper_words


def maybe_naming_repository(shell):
    """The shell once a variable that names the repository may or may not
    be set in it, or set and not exported, set alone (set_alone) or behind
    a time the guard cannot tell from "time" (behind_time_alone): each push
    after it is judged where the shell
    stands, as though nothing were set, and by its words and in the
    session's own clone, as though it were (push_places), so either way the
    clone git runs in is judged."""
    return shell._replace(maybe_named=True)


def moved_to(argv, shell, is_directory):
    """Where the shell stands after the simple command argv (moved_by), in
    each way bash may build its words (word_ways, F29): where the ways
    leave it in different places, as $NOPE cd /x, whose words as written
    run no cd, a shell moved to a directory this module cannot tell, which
    names the repository where every way names it and may where any does."""
    after = [moved_by(way, shell, is_directory) for way in word_ways(argv)[0]]
    if all(each == after[0] for each in after[1:]):
        return after[0]
    named = [each.named_repository for each in after]
    maybe = any(each.maybe_named for each in after) or (any(named) and not all(named))
    return Shell(None, None, True, all(named), maybe)


def moved_by(argv, shell, is_directory):
    """Where the shell stands after the simple command argv: cd, pushd and
    popd, run by the shell itself, behind assignments, reserved words,
    command, time or builtin (IN_THE_SHELL), move it, cd and pushd only to
    a directory (reached), and leave it moved whether they move it or not
    (after_move); GIT_DIR, or another variable that names the repository,
    exported, names a repository this module does not follow
    (naming_repository), and set alone, or behind a time the guard cannot
    tell, may name one (maybe_naming_repository; F29); anything else, a cd another program runs among
    them, as env cd does, leaves it where it was."""
    found = launch(argv)
    if found.start is None:
        named = sets_the_shells_own(found) and any(assigns_repository(word) for word in found.assignments)
        if not named:
            return shell
        if set_alone(found) or behind_time_alone(found):
            return maybe_naming_repository(shell)
        return naming_repository(shell)
    if any(name not in IN_THE_SHELL for name, _ in found.wrappers):
        return shell
    argv = found.argv
    name, operands = argv[found.start], argv[found.start + 1:]
    if name == "builtin" and operands:
        name, operands = operands[0], operands[1:]
    if name in ("export", "declare", "typeset", "readonly", "local"):
        return naming_repository(shell) if any(assigns_repository(word) for word in operands) else shell
    if name == "popd":
        return after_move(shell, popd(shell, operands))
    moves = {"cd": cd, "pushd": pushd}.get(name)
    return after_move(shell, moves(shell, operands, is_directory)) if moves else shell


def placed_commands(command, here, is_directory=os.path.isdir):
    """Each simple command in a shell command, with the shell that runs it
    (Shell), from here, where the command starts: those of the whole
    command, then those of each other reading of it (readings) that no
    reading before placed in the same shell (F25). A move (moved_to) holds
    for what follows it in the same shell and the same reading, a
    subshell's ends at its ), and each part of a pipeline, and a list run
    with & in the background, runs in a shell of its own, so its moves
    change nothing after it. A cd or pushd to what is_directory does not
    take for a directory, by default what is not one as the hook runs,
    leaves the shell where it stood. A directory this module cannot tell
    is None. A reading whose text another before it has is read once,
    since it places the same commands (F29's second review: the reading
    of extended patterns is often the whole reading again)."""
    return placed_in(command, Shell(here), is_directory, 0)


def placed_in(command, start, is_directory, depth):
    """placed_commands's reading of a command from the shell start, at
    depth handings deep: the simple commands of each reading, then those of
    each text the command hands to a shell to run (handed_texts; F32), each
    read once, from the shell that runs it, one level deeper; a handing
    nested past MAX_NESTING is refused unread (HANDED_TOO_DEEP)."""
    placed, handed, known = [], [], set()

    def place(found):
        new = [(argv, shell) for argv, shell in found if (tuple(argv), shell) not in known]
        placed.extend(new)
        known.update((tuple(argv), shell) for argv, shell in new)

    for text in dict.fromkeys(readings(command)):
        place(placed_parts(parts_of(text), start, is_directory, handed))
    for text, shell in dict.fromkeys(handed):
        if not may_push(text):
            continue
        if depth >= MAX_NESTING or too_many_handed(text, shell):
            raise TooDeep(HANDED_TOO_DEEP)
        place(placed_in(text, shell, is_directory, depth + 1))
    return placed


def too_many_handed(text, shell):
    """Whether the text handed to a shell, counted once with the others of
    the run under way (HANDED), passes MOST_HANDED texts, or HANDED_BUDGET
    characters in all (F32's first review: a long text handed on fifteen
    times was read fifteen times and took the guard past its deadline)."""
    seen = HANDED.get()
    if seen is None:
        return False
    if (text, shell) not in seen:
        seen[(text, shell)] = len(text)
        seen[None] = seen.get(None, 0) + len(text)
    return len(seen) - 1 > MOST_HANDED or seen[None] > HANDED_BUDGET


def placed_parts(parts, start, is_directory, handed=None):
    """Each simple command among the parts of one reading of a command
    (parts_of), with the shell that runs it, from the shell start
    (placed_commands); the shell after each simple command found once for
    its words and the shell before it. Each text a simple command hands to
    a shell to run is added to handed, with the shell it runs from
    (handed_texts): what a <<< here-string feeds it, and, at the far end of
    a pipe, what the command before the pipe prints (printed_text)."""
    placed, shell, opened, piped, moves = [], start, [], False, {}
    list_start = pipe_start = shell
    handed = [] if handed is None else handed
    # A reading with no pipe feeds nothing but its here-strings, so its
    # flow is not followed (Flow), which a wide one would pay for.
    flow = Flow() if any(isinstance(part, str) and "|" in part for part in parts) else None
    texts = {}
    for count, part in enumerate(parts):
        if not count % CHECKS_EVERY:
            in_time()
        if isinstance(part, list):
            if flow:
                flow.opens_or_closes(part, before=True)
            fed = standard_input(part, flow.feeding() if flow else ())
            for argv in env_readings(part):
                placed.append((argv, shell))
                handed += [(text, shell) for text in known(texts, handed_texts, argv, fed)]
            if flow:
                flow.printed(known(texts, printed_texts, part, fed))
                flow.opens_or_closes(part, before=False)
            key = (tuple(part), shell)
            if key not in moves:
                moves[key] = moved_to(part, shell, is_directory)
            shell = moves[key]
            continue
        for operator in operators(part):
            if flow:
                flow.operator(operator)
            if operator in ("|", "|&"):
                shell, piped = pipe_start, True
                continue
            if piped:
                shell, piped = pipe_start, False
            if operator == "(":
                opened.append(shell)
            elif operator == ")" and opened:
                shell = opened.pop()
            elif operator == "&":
                shell = list_start
            if operator not in ("&&", "||"):
                list_start = shell
            pipe_start = shell
    return placed


def known(texts, find, argv, fed):
    """What find gives for argv and fed, found once for each pair in one
    reading (placed_parts): a reading of a wide command repeats its simple
    commands by the hundred thousand, and each finding reads its words."""
    key = (find, tuple(argv), fed)
    if key not in texts:
        texts[key] = find(argv, fed)
    return texts[key]


# The words that open a compound command, and those that close one, at
# the head of a simple command's part (Flow), and the reserved words that
# may stand before one there, as then and ! do in then if, ! { and time {
# (F32's second review: only the first word was read, so the second { of
# { { echo ...; }; } | bash opened nothing and the echo's output was lost).
COMPOUND_OPENS = frozenset({"{", "if", "while", "until", "for", "case", "select"})
COMPOUND_CLOSES = frozenset({"}", "fi", "done", "esac"})
COMPOUND_LEADS = frozenset({"!", "time", "then", "do", "else", "elif"})


class Flow:
    """What reaches the standard input of the commands of one reading of a
    command (placed_parts; F32's first review): what the stage before a
    pipe prints, a group's or a subshell's every command's output joined,
    fed to every command of the stage after it, at whatever depth, up to
    the list operator that ends that stage at its own depth."""

    def __init__(self):
        self.depth, self.groups, self.kinds, self.last = 0, [], [], ()
        self.feed, self.feed_depth, self.piped = (), None, False

    def feeding(self):
        return self.feed

    def printed(self, texts):
        self.last = tuple(texts)
        if self.groups:
            self.groups[-1].extend(texts)

    def opens_or_closes(self, part, before):
        self.piped = False
        for word in part:
            if before and word in COMPOUND_OPENS:
                self.open(word)
            elif not before and word in COMPOUND_CLOSES:
                self.close()
            elif word not in COMPOUND_LEADS and word not in COMPOUND_OPENS and word not in COMPOUND_CLOSES:
                break

    def open(self, kind):
        self.depth += 1
        self.groups.append([])
        self.kinds.append(kind)

    def close(self):
        if not self.groups:
            return
        self.depth -= 1
        self.kinds.pop()
        texts = self.groups.pop()
        self.last = ("".join(texts),) if texts else ()
        if self.groups:
            self.groups[-1].extend(self.last)

    def operator(self, operator):
        if self.piped and operator == ";":
            # A newline after a pipe, which the readings write as ;, goes
            # on to the stage after it, as bash reads echo ... |, a newline
            # and bash (F32's second review).
            return
        self.piped = operator in ("|", "|&")
        if self.piped:
            self.feed, self.feed_depth = self.last, self.depth
        elif operator == "(":
            self.open("(")
        elif operator == ")":
            if self.kinds[-1:] != ["case"]:
                # A case's pattern closes with ), which closes no group
                # (F32's second review: case a in a) echo ...;; esac | bash).
                self.close()
        elif self.feed_depth is not None and self.depth <= self.feed_depth:
            self.feed, self.feed_depth = (), None
        if operator not in ("|", "|&", ")"):
            self.last = ()


# The shells a command may hand text to as commands, by name: a -c
# string, or what they read on standard input (handed_texts; F32).
SHELLS = frozenset({"sh", "bash", "rbash", "dash", "ash", "zsh", "ksh", "ksh93", "mksh", "pdksh", "posh", "yash"})
# The options of those shells that take the next word as their value, and
# the long ones; every other option is a flag.
SHELL_VALUE_LETTERS = "oO"
SHELL_VALUE_OPTIONS = frozenset({"--rcfile", "--init-file"})
# The files a shell or source reads its standard input by.
STANDARD_INPUT = frozenset({"/dev/stdin", "/dev/fd/0", "/proc/self/fd/0"})
# The builtins that run their words as commands in the shell: eval its
# words joined, trap its action when a signal comes; and those that run a
# file, read here only where it is standard input.
EVALS = frozenset({"eval"})
TRAPS = frozenset({"trap"})
SOURCES = frozenset({"source", "."})
# The most texts one command may hand to shells to run, all levels
# counted, before the guard refuses it unread: each is read whole, and a
# command near MAX_COMMAND of distinct handings would take the guard past
# its deadline (F32).
MOST_HANDED = 256
HANDED_BUDGET = 2 * MAX_COMMAND
HANDED_TOO_DEEP = (f"a command that hands text to a shell to run (bash -c, eval, trap, a here-string or a pipe "
                   f"into a shell) more than {MAX_NESTING} deep, or more than {MOST_HANDED} texts or "
                   f"{HANDED_BUDGET} characters in all, is refused: the guard cannot read it in time to tell which "
                   "of its commands run; write the text to a file first (.claude/rules/board.md)")
# The texts handed to shells in the guard's run under way (placed_in),
# None outside a run.
HANDED = contextvars.ContextVar("HANDED", default=None)


def handed_texts(argv, fed=()):
    """The texts one reading of a simple command (env_readings) hands to a
    shell to run as commands, each a command of its own (F32): a shell's -c
    string, as in bash -c 'git push origin main'; the words eval joins, and
    trap's action; and, for a shell that reads its standard input, or
    source or . of /dev/stdin, each of fed, the texts that reach its
    standard input, a here-string's or what the command before a pipe
    prints (placed_parts). Until F32 each was one word to the guard, so the
    push in it passed, though a push in a here-document a shell reads was
    refused, since the whole reading reads every body as commands."""
    name, rest = command_words(argv)
    if name in EVALS:
        words = rest[1:] if rest[:1] == ["--"] else rest
        return [" ".join(words)] if words else []
    if name in TRAPS:
        return trap_action(rest)
    if name in SOURCES:
        return list(fed) if rest[:1] and rest[0] in STANDARD_INPUT else []
    if name in SHELLS:
        return shell_texts(rest, fed)
    return []


def command_words(argv):
    """The name of the command a simple command runs, past its assignments
    and wrappers (launch), and the words after it; None and no words where
    it runs none."""
    found = launch(argv)
    if found.start is None or found.start >= len(found.argv):
        return None, []
    words = list(found.argv[found.start:])
    while words and program(words[0]) in BUILTIN:
        words = words[2:] if words[1:2] == ["--"] else words[1:]
    if not words:
        return None, []
    return program(words[0]), words[1:]


# The builtin that runs the builtin its next word names, as builtin eval
# does (F32's first review).
BUILTIN = frozenset({"builtin"})


def trap_action(rest):
    """The action trap sets for its signals, where the words after its
    options set one: trap ACTION SIGNAL..., not trap - SIGNAL..., which
    resets them, nor trap SIGNAL, which resets the one."""
    while rest and rest[0].startswith("-") and rest[0] != "-":
        if rest[0] == "--":
            rest = rest[1:]
            break
        rest = rest[1:]
    return [rest[0]] if len(rest) > 1 and rest[0] != "-" else []


def shell_texts(rest, fed):
    """What a shell given the words rest runs as commands: with -c, its
    first word that is no option; else, where it names no script, or names
    /dev/stdin, or reads its standard input by -s, each of fed. A word of -
    or -- ends the options, as bash reads them; -o and -O take the next
    word, as --rcfile and --init-file do."""
    command, stdin, i = False, False, 0
    while i < len(rest):
        word = rest[i]
        if word in ("-", "--"):
            i += 1
            break
        if word in SHELL_VALUE_OPTIONS:
            i += 2
            continue
        if word.startswith("--") or not word[:1] in ("-", "+") or len(word) < 2:
            if word.startswith("--"):
                i += 1
                continue
            break
        letters = word[1:]
        command = command or (word[0] == "-" and takes_command_string(letters))
        stdin = stdin or (word[0] == "-" and "s" in letters)
        i += 1 + sum(1 for letter in letters if letter in SHELL_VALUE_LETTERS)
    operands = rest[i:]
    if command:
        return operands[:1]
    return shell_input(operands, stdin, fed)


def takes_command_string(letters):
    """Whether a shell's option letters hold -c, which makes its first
    operand the text it runs as commands."""
    return "c" in letters


def shell_input(operands, stdin, fed):
    """What a shell with no -c runs: each of fed, the texts on its
    standard input, where it names no script, names /dev/stdin, or reads
    standard input by -s; nothing where it runs a script, whose text this
    module does not read."""
    if stdin or not operands or operands[0] in STANDARD_INPUT:
        return list(fed)
    return []


def printed_texts(argv, fed=()):
    """What a simple command prints, where this module can tell, for a
    shell at the far end of a pipe (placed_parts): echo's words joined by
    blanks, as written and with its backslash escapes read as $'...' reads
    them (ansi_c_text), since bash's echo reads them only with -e and sh's
    always; printf's format with each %s or %b given the next word, its
    escapes read so, as printf always reads them; and what cat with no file prints,
    its standard input (fed). Nothing for any other command, whose output
    this module does not know."""
    return tuple(text.replace("\0", "") for text in printed_with_nuls(argv, fed))


def printed_with_nuls(argv, fed):
    """What printed_texts gives, a NUL byte still in it, which a shell
    reading the text from a pipe drops (F32's second review: printf
    'gi\\0t push origin main' | bash runs the push). Each text ends as the
    command ends its output: echo with a newline unless -n or \\c ends it,
    printf with none, so a group's outputs join as bash prints them (Flow)."""
    name, rest = command_words(argv)
    if name == "cat" and all(word.startswith("-") for word in rest):
        return tuple(fed)
    if name == "echo":
        ends = "\n"
        while rest and ECHO_OPTIONS.fullmatch(rest[0]):
            ends = "" if "n" in rest[0] else ends
            rest = rest[1:]
        text = " ".join(rest)
        decoded = echo_text(text)
        decoded += "" if ECHO_STOP in text else ends
        return (text + ends,) if decoded == text + ends else (text + ends, decoded)
    if name == "printf":
        if rest[:1] == ["--"]:
            rest = rest[1:]
        if not rest or rest[0].startswith("-v"):
            return ()
        return (printf_text(rest[0], rest[1:]),)
    return ()


def standard_input(part, piped):
    """The texts that reach a simple command's standard input, where this
    module can tell: what its <<< here-string gives (Words.fed), else
    piped, what the command before a pipe prints (printed_texts)."""
    fed = getattr(part, "fed", None)
    return (fed,) if fed is not None else tuple(piped)


ECHO_OPTIONS = re.compile(r"-[neE]+")
# echo -e's escapes, and %b's: \0 with up to three octal digits, where
# printf's format takes \ and one to three; \c ends the output (F32's
# first review).
ECHO_OCTAL = re.compile(r"\\0([0-7]{0,3})")
ECHO_STOP = "\\c"
# A printf conversion: its flags, width and precision, and its letter.
# Its length modifiers, which bash reads past (F32's second review: %ls
# prints as %s), q not among them, since %q is a conversion of its own.
PRINTF_CONVERSION = re.compile(r"%([-+ #0]*)(\*|[0-9]*)(?:\.(\*|[0-9]*))?[hlLjzt]*([a-zA-Z%])")
# The most digits of a width or precision read as a number: one longer is
# read as MAX_COMMAND, past which nothing printed is read (printf_text).
MOST_WIDTH_DIGITS = 9
# The most times printf's format is read again for its words, past which,
# or past MAX_COMMAND characters printed, the rest is not read.
MOST_PRINTF_PASSES = 256


def echo_text(text):
    """What echo -e, and sh's echo always, prints for text: its escapes
    read as $'...' reads them, but \\0nnn as an octal byte, and nothing
    after \\c."""
    text = text.split(ECHO_STOP, 1)[0]
    return escapes_read(ECHO_OCTAL.sub(lambda m: "\\" + (m.group(1) or "0").rjust(3, "0"), text))


def printf_text(form, words):
    """What printf prints for its format and words: the format read again
    while words remain, each conversion given the next word, %s and %b by
    their width and precision, %b's escapes as echo's, %c its first
    character, %q the word quoted, a numeric one the word as written, %%
    a %, and the format's own escapes read as $'...' reads them."""
    out, given, decoded = [], list(words), escapes_read(form)
    for _ in range(MOST_PRINTF_PASSES):
        before = len(given)
        if printf_pass(decoded, given, out) or not given or len(given) == before or printf_size(out) > MAX_COMMAND:
            break
    return "".join(out)


def printf_pass(form, given, out):
    """Read the decoded format once into out, each conversion filled from
    given; whether a %b's \\c stopped printf, which prints nothing more
    (F32's second review: printf '%b%s' 'git push origin main\\c' x | bash
    runs the push). Past MAX_COMMAND characters each padding is cut to one
    blank, so no width builds more than the guard reads (F32's second
    review: printf '%999999999999999s' x took the hook down)."""
    at = 0
    for match in PRINTF_CONVERSION.finditer(form):
        out.append(form[at:match.start()])
        word, stopped = printf_conversion(match, given, MAX_COMMAND - printf_size(out))
        out.append(word)
        if stopped:
            return True
        at = match.end()
    out.append(form[at:])
    return False


def printf_size(out):
    return sum(len(piece) for piece in out)


def printf_conversion(match, given, room=MAX_COMMAND):
    """One printf conversion filled from given, the words left, which it
    takes from, padded to its width or to room, if less, but by one blank
    at least; and whether it is a %b whose \\c stops printf."""
    flags, width, precision, letter = match.groups()
    if letter == "%":
        return "%", False
    for star in (width, precision):
        if star == "*" and given:
            given.pop(0)
    word, stopped = given.pop(0) if given else "", False
    if letter == "b":
        stopped = ECHO_STOP in word
        word = echo_text(word)
    elif letter == "c":
        word = word[:1]
    elif letter == "q":
        word = shlex.quote(word)
    if letter in "bs" and precision and precision.isdigit():
        word = word[:printf_number(precision)]
    if width.isdigit():
        size = min(printf_number(width), max(room, len(word) + 1))
        word = word.ljust(size) if "-" in flags else word.rjust(size)
    return word, stopped


def printf_number(digits):
    """A width or precision as a number, one of more than
    MOST_WIDTH_DIGITS digits as MAX_COMMAND, so int() never meets one too
    long to read (the Constraint on 4,300 digits)."""
    return MAX_COMMAND if len(digits) > MOST_WIDTH_DIGITS else int(digits)


# What a refusal adds where the reading that refused the push is one the
# command's words do not name, which says why it was judged there and how
# to have it judged in its clone alone (push_places); the refusal itself
# reads as it read before F21's fix. UNMOVED is the reading from where the
# command starts, after a cd, pushd or popd, as if none ran, and, where a
# relative -C names the clone, in the session's own clone (after_a_move);
# UNTOLD the session's own clone, where this module cannot tell the clone
# or reads it from a ~ (from_a_tilde).
UNMOVED = ("; the push was judged from where the command starts too, since the hook cannot tell where its cd, pushd "
           "or popd leave the shell: name the clone with git -C and an absolute path that holds nothing the "
           "shell expands, or put no cd, pushd or popd before the push, to have it judged there alone")
UNTOLD = ("; the push was judged in the session's own clone, since the hook cannot tell which clone it runs "
          "in: name the clone with git -C and an absolute path that holds nothing the shell expands, with "
          "nothing naming its repository as GIT_DIR, --git-dir and --bare do, to have the push judged there "
          "alone")


class Reading(NamedTuple):
    """One clone a push may run in, None for its words alone, and what a
    refusal there adds."""

    place: object
    added: str = ""


def followed_and_unmoved(argv, shell, here):
    """Where the git push argv runs as this module follows the shell's
    moves (shell), and where it would run had none of them run, from here
    (push_place); both None, its words alone, while a variable set alone or
    exported before it names its repository (Shell.named_repository),
    which it names whatever cd or -C comes after (F21's fourth review:
    after an export of GIT_DIR an absolute -C was taken for the clone the
    push runs in, alone, though git took the refs GIT_DIR named, so a lease
    onto another worker's claim that 0.5.0's guard refused passed)."""
    if shell.named_repository:
        return None, None
    return push_place(argv, shell.directory), push_place(argv, here)


def after_a_move(argv, shell):
    """Whether the git push argv may run in the session's own clone though
    this module reads it elsewhere: after any cd, pushd or popd, made or
    failed (Shell.moved), where the clone it runs in hangs on where the
    shell stands, as a relative -C does, which is read from where this
    module puts the shell, a cd into a directory the command makes first
    counting as failed and one behind a test that fails as made; an
    absolute -C names its clone wherever the shell stands (F21's fourth
    review: "mkdir -p tmp && cd tmp && git -C .. push origin HEAD" from a
    clone on main was judged in its parent alone, and passed)."""
    return shell.moved and push_place(argv, None) is None


def from_a_tilde(argv):
    """Whether the git push argv names its clone from a ~: of the paths it
    changes directory by (chdir_paths), the last that names a place
    wherever the shell stands, an absolute one or one that starts with a
    ~, starts with a ~, and any after it are joined on what the ~ names.
    This module reads a ~ from the hook's own HOME, while the shell reads
    it from a HOME the command may set before it, and leaves one in
    quotes, which this module does not keep, as written, so such a push
    may run elsewhere than this module reads it, in the session's own
    clone among other places (F21's fourth review: from b-own, which had
    fetched a's claim, "HOME=<the clones>; git -C ~/b-own push
    --force-with-lease origin HEAD:claude/7" was judged in the hook's own
    ~/b-own alone and passed, where 0.5.0's guard refused it)."""
    anchors = [path for path in chdir_paths(argv) or () if path.startswith("~") or os.path.isabs(path)]
    return bool(anchors) and anchors[-1].startswith("~")


def push_places(argv, shell, here):
    """The clones the git push argv may run in, each to judge it in, from
    shell, where the shell stands as this module follows its moves
    (placed_commands), and here, where the command starts: the clone the
    moves take it to, its words alone where this module cannot tell which
    (followed_and_unmoved); then, where it differs, the clone it runs in if
    none of the moves ran, since this module reads a cd in a
    here-document's body, after a test that fails or in a function never
    called as one the shell makes; and here, the session's own clone, as
    every push was judged before F21's fix, where either cannot be told or
    a ~ names the clone (from_a_tilde), or after a move where a relative -C
    names it (after_a_move); and, where a variable that names the
    repository may be set (Shell.maybe_named), here and its words alone
    too. git's and env's -C are joined on in each. Each reading after the
    first is one the command's words do not name, and
    carries what a refusal there adds (UNMOVED, UNTOLD)."""
    followed, unmoved = followed_and_unmoved(argv, shell, here)
    found = [Reading(followed)]
    if unmoved is not None and unmoved != followed:
        found.append(Reading(unmoved, UNMOVED))
    if here is not None and here not in (followed, unmoved):
        if followed is None or unmoved is None or from_a_tilde(argv) or shell.maybe_named:
            found.append(Reading(here, UNTOLD))
        elif after_a_move(argv, shell):
            found.append(Reading(here, UNMOVED))
    if shell.maybe_named and None not in (reading.place for reading in found):
        found.append(Reading(None, UNTOLD))
    return found


def refused_push(command, here, branches_of=None, reader_of=None, is_directory=None):
    """Why the first push of a shell command that this module refuses is
    refused, or None when it refuses none. Each push is judged in each
    clone it may run in, from here, the directory the command starts in
    (placed_commands, push_places), its current and default branch read
    there (branches_of, by default branches) and, for a lease, its claim
    files (reader_of, by default claim_file_reader), each place asked once,
    and refused when any judgement refuses it, with what that reading adds
    (F21: until the fix every push was judged in here's clone, so "git -C
    ../clone push" from a claim branch sent ../clone's main, and a lease
    from a clone that had fetched another's claim was judged in one that
    had not; since its third review a push after a move is judged where
    the command starts too, and one whose clone this module cannot tell in
    here's; and since its fourth one by a relative -C after a move, or by a
    -C from a ~, in here's too, and one whose repository a variable names
    by its words and in here's whatever cd or -C comes after. So every push
    is judged in here's clone, as before the fix, save one whose clone
    git's or env's -C names by an absolute path, or a relative one with no
    cd, pushd or popd before it, holding nothing the shell expands, and
    whose repository nothing names: that push is judged where the -C names
    alone, which is where it runs save where the shell does what this
    module does not read, as the module's docstring says). A push judged by
    its words alone has no current branch and no claim files, the default
    branch being here's. is_directory, by default os.path.isdir, says which
    directories a cd or pushd can reach (placed_commands). The pushes are
    those any reading of the command finds, the whole command's first
    (readings; F25: until its fix a quote that never closed hid every
    push, and a quote in a here-document's body could hide one after it),
    so a command this module read before that fix is refused as it was; a
    command whose bodies nest deeper than MAX_NESTING, or one longer than
    MAX_COMMAND, is refused unread (TOO_DEEP, TOO_LONG), since a hook that
    ran out of time would pass it, and so is one the guard takes longer
    than DEADLINE_SECONDS over, reading it and judging its pushes, git's
    answers included (TOO_SLOW; F25's fourth review: until then the
    deadline held the reading alone, and thousands of pushes, each into a
    clone of its own, ran the judging past the hook's timeout)."""
    if len(command) > MAX_COMMAND:
        return TOO_LONG
    deadline, launches = DEADLINE.set(clock() + DEADLINE_SECONDS), LAUNCHES.set({})
    ways, built, handed = BUILT_WAYS.set({}), READINGS.set({}), HANDED.set({None: 0})
    try:
        placed = placed_commands(command, here, is_directory or os.path.isdir)
        return first_refusal(placed, here, branches_of or branches, reader_of or claim_file_reader)
    except TooDeep as problem:
        return problem.args[0] if problem.args else TOO_DEEP
    except TooSlow:
        return TOO_SLOW
    except (MemoryError, OverflowError, RecursionError, ValueError):
        # A command the guard cannot read for its size refuses, since a hook
        # that fails lets the command run (F32's second review).
        return UNREADABLE
    finally:
        DEADLINE.reset(deadline)
        LAUNCHES.reset(launches)
        BUILT_WAYS.reset(ways)
        READINGS.reset(built)
        HANDED.reset(handed)


def first_refusal(placed, here, branches_of, reader_of):
    """Why the first push among placed, the simple commands of a command
    each with its shell (placed_commands), that a judgement refuses is
    refused, judged in each clone it may run in (push_places), or None; each
    clone is asked once, and only while the run is in time (in_time), and
    each simple command's words are read for a push once, however many
    times the readings place them (F29's timing: a body of 29,000 $(a)$(b)
    read as commands placed 116,000 of them)."""
    asked, pushes = {}, {}

    def shown(place):
        if place not in asked:
            in_time()
            asked[place] = (branches_of(place), reader_of(place))
        return asked[place]

    for argv, shell in placed:
        key = tuple(argv)
        if key not in pushes:
            pushes[key] = push_arguments(argv)
        args = pushes[key]
        if args is None:
            continue
        for place, added in push_places(argv, shell, here):
            (current, default), reader = shown(place) if place is not None else ((None, shown(here)[0][1]), None)
            reason = refusal(args, current, default, reader)
            if reason:
                return reason + added
    return None


def long_option(word):
    """The git push option a word that starts with -- names, as (name,
    negated), or None when git would refuse the word, as an option it does
    not know or a start of more than one, and push nothing."""
    name = word[2:].split("=", 1)[0]
    if name in PUSH_OPTIONS:
        return name, False
    if name.startswith("no-") and name[3:] in PUSH_OPTIONS and name[3:] not in NOT_NEGATED:
        return name[3:], True
    found = {(option, False) for option in PUSH_OPTIONS if option.startswith(name)}
    found |= {(option, True) for option in PUSH_OPTIONS
              if option not in NOT_NEGATED and f"no-{option}".startswith(name)}
    return found.pop() if len(found) == 1 else None


def short_options(bundle, options):
    """Put in options what a bundle of short ones names, -fu being -f and
    -u; 1 when the last one takes the next word as its value, as -o does,
    else 0. A letter git does not know makes it refuse the push."""
    for at, letter in enumerate(bundle):
        name = PUSH_LETTERS.get(letter)
        if name is None:
            return 0
        options.pop(name, None)
        options[name] = True
        if PUSH_OPTIONS[name] == VALUE:
            return 1 if at == len(bundle) - 1 else 0
    return 0


def push_words(args):
    """git push's arguments read as git reads them: the options in force,
    by name, in the order they were last given, and the positional words,
    the remote and then the refspecs. --force-with-lease is in force as
    the list of what each giving of it named, "" for one given bare,
    which holds every ref pushed to its remote-tracking ref, and
    "ref" or "ref:commit" for one that names a ref (git push's --help);
    --no-force-with-lease takes them all back, as git does."""
    options, positionals, i = {}, [], 0
    while i < len(args):
        word = args[i]
        i += 1
        if word == "--":
            positionals += args[i:]
            break
        if word.startswith("--"):
            found = long_option(word)
            if found is None:
                continue
            name, negated = found
            named = options.pop(name, None)
            if negated:
                continue
            if PUSH_OPTIONS[name] == OPTIONAL:
                options[name] = (named if isinstance(named, list) else []) + [word.partition("=")[2]]
            else:
                options[name] = True
                if PUSH_OPTIONS[name] == VALUE and "=" not in word:
                    i += 1
        elif word.startswith("-") and word != "-":
            i += short_options(word[1:], options)
        else:
            positionals.append(word)
    return options, positionals


def moves_many(spec):
    """Whether a refspec pushes many branches: the matching refspec :,
    which updates every branch both sides have, or a pattern with a *."""
    return spec.lstrip("+") == ":" or "*" in spec


def branch_named(name):
    """The branch one side of a refspec names, as git resolves it: a name
    with refs/ left off is read as one under refs/, so heads/main is
    refs/heads/main, which is main, on either side of the colon (git 2.43
    pushed main on "git push origin heads/main" and on "HEAD:heads/main"),
    after the colon where the remote has the branch, as it always has the
    default branch and board, git making heads/X of a name it lacks.
    refs/heads/ comes off once, so refs/heads/heads/main is the branch
    heads/main; any other name under refs/, a tag's, is kept whole and is
    no branch the rules name."""
    if name.startswith("refs/"):
        return name.removeprefix("refs/heads/")
    return name.removeprefix("heads/")


def refspec_pairs(refspecs, current_branch, deleting):
    """What a push writes, as (source, branch, deleted) for each refspec: the
    branch is the refspec's destination, after the colon when it has one,
    as git resolves the name (branch_named), the current branch for HEAD,
    for @ alone, which git reads as HEAD, or for no refspec at all; the
    source is what stands before the colon, else the refspec itself, HEAD
    for @ alone and for no refspec, as git names what it would push; and
    a deletion is one under --delete or a refspec with nothing before its
    colon."""
    pairs = []
    for spec in refspecs:
        name, deleted = spec.lstrip("+"), deleting
        source = name
        if ":" in name:
            source, name = name.split(":", 1)
            deleted = deleted or not source
        elif name == "@":
            source = name = "HEAD"
        name = branch_named(name)
        if name in ("", "HEAD"):
            name = current_branch
        if name:
            pairs.append((source, name, deleted))
    if not refspecs and current_branch:
        pairs.append(("HEAD", current_branch, deleting))
    return pairs


def targets_of(refspecs, current_branch, deleting):
    """The branches a push writes, each with whether it deletes it
    (refspec_pairs)."""
    return [(name, deleted) for _, name, deleted in refspec_pairs(refspecs, current_branch, deleting)]


def misnamed_claim(name):
    """Why a push that writes this branch is refused for a name that reads
    as a claim and claims nothing (CLAIM_NUMBER), naming the issue's own
    claim branch where the digits spell one; None for any other name."""
    key = name[len("claude/"):] if name.startswith("claude/") else ""
    if not key.isdigit() or CLAIM_NUMBER.fullmatch(key):
        return None
    meant = "".join(str(unicodedata.digit(digit)) for digit in key).lstrip("0") or "0"
    own = f", so the claim branch of #{meant} is claude/{meant}" if len(meant) <= 9 else ""
    return (f"a push to {name} is refused: it reads as a claim branch and claims nothing, since a claim "
            f"branch is claude/ and the issue's number in one to nine ASCII digits with no leading "
            f"zero{own} (docs/agents/PROTOCOL.md, Claim)")


def claim_number(name):
    """The issue a branch claims, as text, or None for any other name: the
    protocol's Claim rule, CLAIM_NUMBER after claude/."""
    key = name[len("claude/"):] if name.startswith("claude/") else ""
    return key if CLAIM_NUMBER.fullmatch(key) else None


def claim_path(number):
    return f"{CLAIM_DIR}/{number}.md"


def worker_named(text):
    """The worker a claim file names (WORKER_LINE), or None when no line of
    it does, or there is no file."""
    for line in (text or "").splitlines():
        match = WORKER_LINE.match(line)
        if match:
            return TRAILING_NOTE.sub("", match.group(1)).strip("`*_ ")
    return None


def lease_holds(values, name):
    """How a --force-with-lease given as values (push_words) holds the
    branch name: None when it does not, so git pushes that ref as a plain
    push would; else the commit the lease expects there, "" for the
    remote-tracking ref, as a bare --force-with-lease, which holds every
    ref not named, or one that names the ref alone expects it. A giving
    that names the ref holds it whatever a bare one says, the first such
    giving winning, as git applies the first entry that matches the ref."""
    for value in values:
        ref, _, commit = value.partition(":")
        if ref and branch_named(ref) == name:
            return commit
    return "" if "" in values else None


def remote_ref(remote, name):
    """The remote-tracking ref of a branch for the remote a push names,
    a name with no slash or colon; a URL or a path has none, and origin's
    stands in, as the clone's best knowledge of the branch."""
    named = remote if remote and "/" not in remote and ":" not in remote else "origin"
    return f"{REMOTE_REFS}{named}/{name}"


def other_holder(name, source, expected, claim_file_at, remote="origin"):
    """Why a lease onto the claim branch name is refused as another
    worker's: the claim file on the branch the lease holds against
    names a worker, and the claim file on the branch pushed, source,
    names another or is not there. The branch the lease holds against is
    the commit expected, when the lease names one git can show, else the
    remote-tracking ref (remote_ref): a lease written with a shell
    variable names a commit git cannot show, and the push lands only when
    the remote stands at the commit meant, of which the remote-tracking
    ref is the clone's best knowledge. claim_file_at gives a claim file's
    text by a ref and an issue number, or None (claim_file_reader). None
    when the lease is allowed, or when the fence cannot judge it."""
    number = claim_number(name)
    if number is None or expected is None:
        return None
    held = worker_named(claim_file_at(expected, number)) if expected else None
    if held is None:
        held = worker_named(claim_file_at(remote_ref(remote, name), number))
    if held is None:
        return None
    pushed = worker_named(claim_file_at(source, number))
    if pushed == held:
        # The holder's own rebase, or a branch built on the holder's own
        # claim commit, which keeps the holder's work and name: the fence
        # judges whose claim a branch carries, not who pushes, since the
        # session's own name is nowhere in the command (G11b).
        return None
    theirs = f"{name} on origin is {held}'s claim (its claim file {claim_path(number)} says so)"
    yours = f"claims #{number} for {pushed}" if pushed else f"carries no claim file for #{number}"
    return (f"a push with --force-with-lease to {name} is refused: {theirs} and the branch you push {yours}; "
            f"never push to another worker's branch (.claude/rules/board.md)")


# The refusals of a push whose option or branch bash builds (F29): an
# option whose name bash builds may be --force or -f, and a branch whose
# name it builds around its own characters, as ma$(true)in or mai?, may be
# the default branch or board. A branch that is one expansion whole, as
# "$(git branch --show-current)", is a value the guard cannot know, read
# as a variable's is, and so is one made of expansions alone, as
# "${PREFIX}${N}" (expansions_alone), and as nothing where bash may make it
# nothing (unset_value).
BUILT_OPTION = ("a push with the option {word} is refused: the shell builds its name, which may be --force; write "
                "the option out (.claude/rules/board.md)")
BUILT_TARGET = ("a push to {word} is refused: the shell builds that name, which may be {meant}; write the branch "
                "out (.claude/rules/board.md)")


def as_shown(word):
    """A word as a refusal names it, each substitution read apart shown as
    $(...) where the reading holds SUBSTITUTED, a character no terminal
    shows, so the refusal stays ASCII as the words written are."""
    return word.replace(SUBSTITUTED, "$(...)")


def built_option(args):
    """The first of a push's option words whose name, before any =, bash
    builds (is_built), as -$(echo f), --for${X}ce or --"$opt", else None,
    unless it opens with a whole option no other extends (whole_flag); and
    the first word that does not open with - as written but that bash may
    build into a forcing option (may_force), as $(echo -)f or a * a file
    named -f fills (F29's second review); words after a -- are no
    options."""
    for word in args:
        if word == "--":
            return None
        if word.startswith("-") and is_built(word.partition("=")[0]) and not whole_flag(word):
            return word
        if not word.startswith("-") and may_force(word):
            return word
    return None


# How may_force reads a word: nothing yet, then a - read with no f after
# it, then a short bundle that holds an f.
NOTHING_YET, DASHED, FORCING = 0, 1, 2


def may_force(word):
    """Whether bash may build a word, which is built and not made of
    expansions alone (expansions_alone), into --force or a short bundle
    that holds f, as -f and -uf are: its own characters, in order, around
    holes filled as they may be (word_pieces, hole_allows). A word of
    expansions alone, as "$BRANCH", is a value the guard cannot know, read
    as before."""
    if not is_built(word) or expansions_alone(word):
        return False
    if may_be(word, "--force"):
        return True
    states = {NOTHING_YET}
    for piece in word_pieces(word):
        if isinstance(piece, str):
            for ch in piece:
                states = {state for state in map(lambda at: bundle_step(at, ch), states) if state is not None}
        elif piece.any:
            states |= {DASHED, FORCING}
        else:
            states = {state for at in states for ch in "-fv" if hole_allows(piece, ch)
                      for state in [bundle_step(at, ch)] if state is not None}
        if not states:
            return False
    return FORCING in states


def bundle_step(state, ch):
    """Where may_force stands after the character ch: a - opens a bundle,
    a letter goes on with one, an f makes it force; None where the word can
    no longer be one."""
    if state == NOTHING_YET:
        return DASHED if ch == "-" else None
    if not ch.isascii() or not ch.isalpha():
        return None
    return FORCING if ch == "f" or state == FORCING else DASHED


def may_open_with_plus(spec):
    """Whether a refspec's first part is a hole of a pattern that may give
    a +, as ?main and [+]main do where a file is named +main (F29's second
    review): a + forces that ref. An expansion's hole is a value the guard
    cannot know, as $B is, read where it is nothing by the unset world."""
    pieces = word_pieces(spec) if is_built(spec) else []
    return bool(pieces) and isinstance(pieces[0], Hole) and not pieces[0].expansion and hole_allows(pieces[0], "+")


def colon_given(spec, branches):
    """What a refspec with no colon of its own may be where a hole of a
    pattern in it gives one, as a file's name may hold one (F29's second
    review): the refspec :, which moves every branch both sides have, as ?
    and [:] are where a file is named :, or one that writes or deletes a
    branch of branches, as [:]main is where a file is named :main; else
    None. One with an expansion in it is read where the expansion is
    nothing by the unset world."""
    if not is_built(spec):
        return None
    pieces = word_pieces(spec)
    if any(isinstance(piece, Hole) and piece.expansion for piece in pieces):
        return None
    if pieces_may_be(pieces, ":"):
        return ":"
    if any(isinstance(piece, str) and ":" in piece for piece in pieces):
        return None
    for at, piece in enumerate(pieces):
        if isinstance(piece, Hole) and hole_allows(piece, ":"):
            meant = next((branch for branch in branches if branch and colon_then(piece, pieces[at + 1:], branch)), None)
            if meant:
                return meant
    return None


def colon_then(hole, after, branch):
    """Whether the hole may give a colon and the pieces after it the name
    of branch, as git resolves it (branch_named): a * may give the start of
    the name with the colon, a ? or a class the colon alone."""
    for name in (branch, "heads/" + branch, "refs/heads/" + branch):
        starts = range(len(name) + 1) if hole.any else (0,)
        if any(pieces_may_be(after, name[start:]) for start in starts):
            return True
    return False


def whole_flag(word):
    """Whether a long option word opens, before the first part bash builds,
    with the whole name of an option of git push that takes no value and
    that no other option's name extends, as --no-verify${NOPE} does (F29's
    first review). An unquoted hole may split into more words, so bash may
    build more than that option; but where the hole's text opens with a
    blank, what it makes is what it makes of the option followed by the
    hole written as a word of its own, --no-verify $NOPE, which already
    passes, an expansion whole being a value the guard cannot know, and
    where it does not, the option's name runs on into a name no option of
    git push has, since none extends it, and git refuses the push (F29's
    second review). A built short bundle is never one, since -u${X} may be
    -uf, and --"$opt" is none, since it may be --force, nor is --force$X,
    which may be --force-with-lease."""
    head = next(iter(word_pieces(word)), "")
    if not isinstance(head, str) or not head.startswith("--"):
        return False
    name = head[2:]
    return name in FLAG_NAMES and not any(other != name and other.startswith(name)
                                          for other in FLAG_NAMES | VALUE_NAMES)


def built_target(name, default_branch):
    """The branch the rules name that bash may build the branch name into,
    the default branch or board (may_be), where the name is built and is
    not made of expansions alone (expansions_alone), else None."""
    if not is_built(name) or expansions_alone(name):
        return None
    return next((branch for branch in (default_branch, "board") if branch and may_be(name, branch)), None)


def refusal(args, current_branch, default_branch, claim_file_at=None):
    """Why this push is refused, or None when it is allowed: each reading
    of Readings, the first refused, and an Unreadable word's own why.
    claim_file_at
    reads a claim file by a ref and an issue number, for the fence on a
    lease onto another worker's claim (other_holder); None, as a caller
    with no clone to read gives, leaves that fence out."""
    if isinstance(args, Readings):
        judged = (refusal(list(each), current_branch, default_branch, claim_file_at)
                  for each in [args, *args.others])
        return next((reason for reason in judged if reason), None)
    unread = next((word for word in args if isinstance(word, Unreadable)), None)
    if unread is not None:
        return str(unread)
    built = built_option(args)
    if built:
        return BUILT_OPTION.format(word=as_shown(built))
    options, positionals = push_words(args)
    remote, refspecs, deleting = (positionals or [""])[0], positionals[1:], "delete" in options
    wholesale = [name for name in options if name in WHOLESALE]
    if wholesale:
        return (f"a push with --{wholesale[-1]} moves branches wholesale; "
                "push your own claim branch alone (.claude/rules/board.md)")
    colons = {} if deleting else {spec: colon_given(spec, (default_branch, "board")) for spec in refspecs}
    many = [spec for spec in refspecs if moves_many(spec) or colons.get(spec) == ":"]
    if many:
        return (f"a push with the refspec {many[0]} moves branches wholesale; "
                "push your own claim branch alone (.claude/rules/board.md)")
    if "force" in options:
        return ("--force overwrites history another account may hold; the one forced push "
                "the protocol allows is --force-with-lease on your own claude/ branch "
                "after a rebase (.claude/rules/board.md)")
    # A refspec's leading + forces that one ref, as --force forces them all,
    # and git lets it override --force-with-lease; the first word is the
    # remote, never a refspec, and under --delete each refspec is a name to
    # delete, the + in it part of the name.
    if not deleting and any(spec.startswith("+") or may_open_with_plus(spec) for spec in refspecs):
        return ("a + before a refspec forces that ref as --force does, and git lets it override "
                "--force-with-lease; the one forced push the protocol allows is --force-with-lease "
                "on your own claude/ branch after a rebase, with no + (.claude/rules/board.md)")
    given = next((spec for spec in refspecs if colons.get(spec)), None)
    if given:
        return BUILT_TARGET.format(word=as_shown(given), meant=colons[given])
    targets = targets_of(refspecs, current_branch, deleting)
    for name, deleted in targets:
        if name == default_branch:
            return (f"a push to {name} is refused: the default branch takes only "
                    "pull requests (.claude/rules/board.md)")
        if name == "board":
            return ("a push to board is refused: only the board workflow writes that "
                    "branch (.claude/rules/board.md)")
        misnamed = None if deleted else misnamed_claim(name)
        if misnamed:
            return misnamed
        meant = built_target(name, default_branch)
        if meant:
            return BUILT_TARGET.format(word=as_shown(name), meant=meant)
    if "force-with-lease" in options:
        strays = [name for name, _ in targets if not name.startswith("claude/")]
        if strays:
            return (f"--force-with-lease to {', '.join(strays)} is refused: a forced push "
                    "is for your own claude/ branch only (.claude/rules/board.md)")
        for source, name, deleted in refspec_pairs(refspecs, current_branch, deleting):
            if deleted or claim_file_at is None:
                continue
            reason = other_holder(name, source, lease_holds(options["force-with-lease"], name), claim_file_at,
                                  remote)
            if reason:
                return reason
    return None


def git_line(args, cwd):
    """One line from git, or None when git cannot answer, as in a directory
    or with an argument that holds a NUL byte or a character the system
    cannot encode, which subprocess refuses with ValueError: the directory
    comes from the command since F21's fix, and the ref of a lease came
    from it before. Inside a guard's run git has what is left of its
    deadline, at most 10 seconds, and a call cut short by the deadline is
    TooSlow, never None, which would let the push pass (DEADLINE_SECONDS)."""
    in_time()
    left = seconds_left()
    timeout = 10 if left is None else min(10, max(left, 0.001))
    try:
        done = subprocess.run(["git", *args], cwd=cwd, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        if timeout < 10:
            raise TooSlow()
        return None
    except (OSError, ValueError):
        return None
    if done.returncode != 0:
        return None
    return done.stdout.decode("utf-8", "replace").strip() or None


def claim_file_reader(cwd):
    """A reader of claim files in the clone at cwd, for the fence on a
    lease (other_holder): the text of docs/agents/claims/N.md at a ref,
    at its tip, else at the commit where the ref's history last added it,
    since a finished branch has deleted the file; None when git shows
    none, for a ref the clone lacks or a branch whose history never had
    the file among other things, which the fence then does not judge.
    The ref comes from the command, so it is handed to git after
    --end-of-options, where one that starts with a dash is a ref and not
    an option."""
    def read(ref, number):
        path = claim_path(number)
        text = git_line(["show", "--end-of-options", f"{ref}:{path}"], cwd)
        if text is None:
            added = git_line(["log", "--format=%H", "--diff-filter=A", "-n", "1", "--end-of-options", ref, "--", path],
                             cwd)
            text = git_line(["show", "--end-of-options", f"{added}:{path}"], cwd) if added else None
        return text
    return read


def branches(cwd):
    """The current branch (None when detached) and the default branch."""
    current = git_line(["rev-parse", "--abbrev-ref", "HEAD"], cwd)
    if current == "HEAD":
        current = None
    head = git_line(["symbolic-ref", "--quiet", "--short", "refs/remotes/origin/HEAD"], cwd)
    default = head.removeprefix("origin/") if head else "main"
    return current, default


def command_of(payload):
    """The shell command a hook payload carries, or "" when the payload is
    not the shape Claude Code sends (an object whose tool_input is an
    object with a string command): a payload the guard cannot read passes,
    since branch protection still stands behind it."""
    tool_input = payload.get("tool_input") if isinstance(payload, dict) else None
    command = tool_input.get("command") if isinstance(tool_input, dict) else None
    return command if isinstance(command, str) else ""


def where_to_ask(payload, cwd):
    """The directory git is asked in: the session's own, from the payload,
    else the caller's, else here. A payload cwd that is not text is passed
    over, since git cannot be run in it."""
    session = payload.get("cwd") if isinstance(payload, dict) else None
    return (session if isinstance(session, str) and session else None) or cwd or "."


def may_push(command):
    """Whether the hook reads the command at all (run): its text, each line
    continuation out, may spell push (holds_push), or bash may build a part
    of a word of it (BUILDS_WORDS), or it may hand a shell a text whose
    escapes or printf's conversions spell push (HANDS_TEXT; F32's first
    review: printf 'git pu\\163h origin main' | bash)."""
    joined = lines_joined(command)
    return (holds_push([joined]) or not BUILDS_WORDS.isdisjoint(command)
            or not BUILDS_TEXT.isdisjoint(command) and bool(HANDS_TEXT.search(joined)))


# A word that may hand a shell a text to run (handed_texts), and what in a
# text it hands lets echo or printf print what is not written out. It is
# looked for once the continuations are out, as bash reads words, and
# never right after one, so the sh that git pu, a continuation and sh
# leaves is no shell's name and the gate reads that push only by joining
# its lines (the mutant "the hook's gate read with continuations in").
HANDS_TEXT = re.compile(r"(?<![\w.-])(?<!\\\n)(?:\w*sh[0-9]*|eval|trap|source|\.)(?![\w.-])")
BUILDS_TEXT = frozenset("\\%")


def lines_joined(command):
    """The command with each line continuation taken out, as bash takes it
    out before it reads a word (run)."""
    return command.replace("\\\n", "")


# What in a command's text lets bash build a word whose text is not written
# out (run): an expansion or a substitution, a backquote, a brace and a
# pattern's characters, an extended pattern's parenthesis among them.
BUILDS_WORDS = frozenset("$`{*?[(")


def run(payload, cwd=None):
    """The hook's answer to one tool call: the exit code and the reason,
    which is empty when the call passes. The payload is the hook's JSON
    already decoded; cli.py reads it from stdin and lets a payload that is
    not JSON pass, as this lets one that is not an object pass. The
    command starts in the payload's own cwd, else in cwd, else here, and
    each push is judged in each clone it may run in from there
    (refused_push), git being asked only once the command holds a push, so
    a session that never pushes never waits on it; and for the claim files
    a lease is judged by, only on a lease push (other_holder). A command
    holds a push where its text spells push once its quotes, backslashes
    and ${NAME}s are taken out (holds_push), as git pu""sh and env -S 'git
    pu""sh ...' do, which bash and env run as push (F28's third review);
    the template's own hook, until it calls this, asks only of a text that
    holds push as written. And it may push wherever bash builds a part of a
    word (BUILDS_WORDS), as git pu$(true)sh, git $'pu\\x73h', git p{u,}sh and
    git pus? do, whose text spells no push (F29), so such a command is read
    too, git still asked only once a push is found in it; and a line
    continuation is taken out first, as bash takes it out, so git p\\, a
    newline and ush is read as the push it runs (F29's first review)."""
    command = command_of(payload)
    if not may_push(command):
        return PASS, ""
    reason = refused_push(command, where_to_ask(payload, cwd))
    return (REFUSE, reason) if reason else (PASS, "")

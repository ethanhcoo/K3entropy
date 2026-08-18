#!/usr/bin/env python3
"""Check that every published copy of the Flipper word agrees."""

import argparse
import re
import subprocess
from pathlib import Path


EXPECTED_GENERATORS = 250
MACHINE_TOKEN = re.compile(r"[sS]_[0-8]\Z")
MANUSCRIPT_LABEL = r"\label{eq:f2-flipper-word}"
MANUSCRIPT_GENERATOR = re.compile(
    r"""
    \\hat\s*(?:\{\s*s\s*\}|s)
    \s*_\s*(?:\{\s*(?P<braced>[0-8])\s*\}|(?P<plain>[0-8]))
    (?P<inverse>\s*\^\s*\{\s*-\s*1\s*\})?
    """,
    re.VERBOSE,
)


def fail(message):
    raise SystemExit(f"mapping-class word check failed: {message}")


def parse_machine_word(data, source):
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError:
        fail(f"{source} is not ASCII")

    if not text.endswith("\n") or text.count("\n") != 1 or "\r" in text:
        fail(f"{source} must be one line terminated by exactly one LF")

    body = text[:-1]
    if not body or body != body.strip() or body.endswith("."):
        fail(f"{source} has leading/trailing whitespace or a trailing dot")

    tokens = body.split(".")
    invalid = [token for token in tokens if MACHINE_TOKEN.fullmatch(token) is None]
    if invalid:
        fail(f"{source} contains invalid token {invalid[0]!r}")
    if len(tokens) != EXPECTED_GENERATORS:
        fail(
            f"{source} has {len(tokens)} generators; "
            f"expected {EXPECTED_GENERATORS}"
        )
    return tokens


def compare(expected, actual, actual_source):
    for index, (left, right) in enumerate(zip(expected, actual), start=1):
        if left != right:
            fail(
                f"{actual_source} first differs at generator {index}: "
                f"expected {left}, found {right}"
            )
    if len(expected) != len(actual):
        fail(
            f"{actual_source} has {len(actual)} generators; "
            f"expected {len(expected)}"
        )


def check_mclass(executable, expected):
    completed = subprocess.run(
        [str(executable), "--flipper-word"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        diagnostic = completed.stderr.decode("utf-8", errors="replace").strip()
        fail(
            f"{executable} --flipper-word exited with status "
            f"{completed.returncode}: {diagnostic}"
        )
    if completed.stderr:
        diagnostic = completed.stderr.decode("utf-8", errors="replace").strip()
        fail(f"{executable} --flipper-word wrote to stderr: {diagnostic}")

    actual = parse_machine_word(completed.stdout, f"{executable} output")
    compare(expected, actual, f"{executable} output")
    print(
        f"PASS: {executable} reproduces the canonical "
        f"{len(expected)}-generator word"
    )


def manuscript_tokens(path):
    text = path.read_text(encoding="utf-8")
    if text.count(MANUSCRIPT_LABEL) != 1:
        fail(
            f"{path} must contain exactly one "
            f"{MANUSCRIPT_LABEL} anchor"
        )

    label_position = text.index(MANUSCRIPT_LABEL)
    begin_position = text.rfind(r"\begin{equation}", 0, label_position + 1)
    end_position = text.find(r"\end{equation}", label_position)
    if begin_position < 0 or end_position < 0:
        fail(f"{MANUSCRIPT_LABEL} is not inside an equation in {path}")

    equation = text[begin_position:end_position]
    tokens = []
    for match in MANUSCRIPT_GENERATOR.finditer(equation):
        index = match.group("braced") or match.group("plain")
        prefix = "S" if match.group("inverse") else "s"
        tokens.append(f"{prefix}_{index}")

    if len(tokens) != EXPECTED_GENERATORS:
        fail(
            f"the equation anchored by {MANUSCRIPT_LABEL} in {path} "
            f"has {len(tokens)} generators; expected {EXPECTED_GENERATORS}"
        )
    return tokens


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--word", required=True, type=Path)
    parser.add_argument("--mclass", type=Path)
    parser.add_argument("--manuscript", type=Path)
    args = parser.parse_args()

    if args.mclass is None and args.manuscript is None:
        parser.error("supply --mclass, --manuscript, or both")

    expected = parse_machine_word(
        args.word.read_bytes(),
        str(args.word),
    )
    if args.mclass is not None:
        check_mclass(args.mclass, expected)
    if args.manuscript is not None:
        actual = manuscript_tokens(args.manuscript)
        compare(expected, actual, str(args.manuscript))
        print(
            f"PASS: {args.manuscript} matches the canonical "
            f"{len(expected)}-generator word"
        )


if __name__ == "__main__":
    main()

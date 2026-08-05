#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT="$DIR/.."

# Mirrors the lint job in .github/workflows/ci.yml. Keep the two in step, or this passes
# locally while CI fails.
#
# The three mypy runs stay separate because petbox/ is a namespace package with no
# __init__.py: one invocation over all three trees makes mypy resolve petbox/dca/__init__.py
# as both `dca` and `petbox.dca` and abort.

echo ruff check petbox/dca tests docs
ruff check "$ROOT/petbox/dca" "$ROOT/tests" "$ROOT/docs"
echo

echo mypy --strict petbox/dca
( cd "$ROOT" && mypy --strict petbox/dca )
echo

echo mypy --strict tests
( cd "$ROOT" && mypy --strict tests )
echo

echo mypy --strict docs
( cd "$ROOT" && mypy --strict docs )
echo

echo pytest --cov=petbox.dca --cov-report=term-missing --hypothesis-show-statistics -v .
pytest --cov=petbox.dca --cov-report=term-missing --hypothesis-show-statistics -v .

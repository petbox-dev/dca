#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
ROOT="$DIR/.."

# Mirrors the lint job in .github/workflows/ci.yml. Keep the two in step, or this passes
# locally while CI fails. mypy runs from the repo root and uses -p test for the test tree:
# test/ is a package, and passing the path makes mypy see test/data.py as both `data`
# and `test.data`, which is a hard error.

echo ruff check petbox/dca test
ruff check "$ROOT/petbox/dca" "$ROOT/test"
echo

echo mypy petbox/dca
( cd "$ROOT" && mypy petbox/dca )
echo

echo mypy --strict petbox/dca
( cd "$ROOT" && mypy --strict petbox/dca )
echo

echo mypy --strict -p test
( cd "$ROOT" && mypy --strict -p test )
echo

echo pytest --cov=petbox.dca --cov-report=term-missing --hypothesis-show-statistics -v .
pytest --cov=petbox.dca --cov-report=term-missing --hypothesis-show-statistics -v .

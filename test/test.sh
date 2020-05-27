#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


echo flake8 ../petbox/dca
flake8 $DIR/../petbox/dca
echo

echo mypy ../petbox/dca
mypy $DIR/../petbox/dca
echo

echo pytest --cov=petbox.dca --cov-report=term-missing --hypothesis-show-statistics -v .
pytest --cov=petbox.dca --cov-report=term-missing --hypothesis-show-statistics -v .

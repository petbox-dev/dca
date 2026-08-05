:: Run tests and generate report
::
:: Mirrors the lint job in .github/workflows/ci.yml. Keep the two in step, or this passes
:: locally while CI fails. mypy runs from the repo root and uses -p test for the test tree:
:: test\ is a package, and passing the path makes mypy see test\data.py as both `data`
:: and `test.data`, which is a hard error.

pushd %~dp0..

ruff check petbox\dca test docs
mypy petbox\dca
mypy --strict petbox\dca
mypy --strict -p test

popd

pytest

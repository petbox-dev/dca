:: Run tests and generate report
::
:: Mirrors the lint job in .github/workflows/ci.yml. Keep the two in step, or this passes
:: locally while CI fails.
::
:: The three mypy runs stay separate because petbox\ is a namespace package with no
:: __init__.py: one invocation over all three trees makes mypy resolve petbox\dca\__init__.py
:: as both `dca` and `petbox.dca` and abort.

pushd %~dp0..

ruff check petbox\dca tests docs
mypy --strict petbox\dca
mypy --strict tests
mypy --strict docs

popd

pytest

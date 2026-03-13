:: Run tests and generate report

ruff check %~dp0..\petbox\dca
mypy %~dp0..\petbox\dca

pytest

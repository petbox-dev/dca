:: Run tests and generate report

flake8 %~dp0..\petbox\dca
mypy %~dp0..\petbox\dca

pytest --cov=%~dp0..\petbox.dca --cov-report=term-missing --hypothesis-show-statistics -v .

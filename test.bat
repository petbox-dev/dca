:: Run tests and generate HTML report

flake8 petbox\dca
mypy petbox\dca

pytest --cov=petbox\dca --cov-report=html -v test
REM coveralls

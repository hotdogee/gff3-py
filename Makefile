.PHONY: clean-pyc clean-build clean-test docs clean lint test test-all coverage dist release

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "lint - check style with flake8"
	@echo "test - run tests quickly with the default Python"
	@echo "test-all - run tests on every Python version with tox"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation, including API docs"
	@echo "dist - build the sdist and wheel and check them with twine"
	@echo "release - upload dist/ to PyPI manually (normally done by the publish workflow)"

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

lint:
	flake8 gff3 tests

test:
	python -m pytest

test-all:
	tox

coverage:
	coverage run --source gff3 -m pytest
	coverage report -m
	coverage html
	open htmlcov/index.html

docs:
	rm -f docs/gff3.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ gff3
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	open docs/_build/html/index.html

dist: clean
	python -m build
	python -m twine check --strict dist/*
	ls -l dist

release: dist
	python -m twine upload dist/*

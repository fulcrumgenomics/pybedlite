# Contributing to pybedlite

## Getting Setup for Development Work

Clone the repository to your local machine.
Note that pybedlite includes [cgranges][cgranges-link] as a submodule, so you must use the `--recurse-submodules` option:

```
git clone --recurse-submodules https://github.com/fulcrumgenomics/pybedlite.git
```

[Poetry][poetry-link] is used to manage the python development environment.

A simple way to create an environment with the desired version of python and poetry is to use [conda][conda-link].  E.g.:

```bash
conda create -n pybedlite python=3.8 poetry=1.6
conda activate pybedlite
poetry install
```

If the methods listed above do not work try the following:

```bash
mamba create -n pybedlite -c conda-forge "python=3.9.16" "poetry=1.6.1"
mamba activate pybedlite
poetry install
```

[poetry-link]: https://github.com/python-poetry/poetry
[conda-link]:  https://docs.conda.io/en/latest/miniconda.html
[cgranges-link]: https://github.com/lh3/cgranges

## Checking the Build

Check the library with:

```bash
poetry run ./ci/check.sh
```

## Building the Documentation

Build the documentation with:

```bash
(cd docs/ && poetry run make html)
```

## Creating a Release on PyPi

1. Clone the repository recursively and ensure you are on the `main` (un-dirty) branch
2. Checkout a new branch to prepare the library for release
3. Bump the version of the library to the desired SemVer with `poetry version #.#.#`
4. Commit the version bump changes with a Git commit message like `chore(release): bump to #.#.#`
5. Push the commit to the upstream remote, open a PR, ensure tests pass, and seek reviews
6. Squash merge the PR
7. Tag the new commit on the main branch of the repository with the new SemVer

GitHub Actions will take care of the remainder of the deploy and release process with:

1. Unit tests will be run for safety-sake
2. A source distribution will be built
3. Many multi-arch multi-Python binary distributions will be built
4. Assets will be deployed to PyPi with the new SemVer
5. A [Conventional Commit](https://www.conventionalcommits.org/en/v1.0.0/)-aware changelog will be drafted
6. A GitHub release will be created with the new SemVer and the drafted changelog

Consider editing the changelog if there are any errors or necessary enhancements.

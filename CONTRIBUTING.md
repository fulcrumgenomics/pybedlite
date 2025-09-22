# Development and Testing

## Development Installation

Install the Python package and dependency management tool [`uv`](https://docs.astral.sh/uv/getting-started/installation/) using official documentation.

Install the dependencies of the project with:

```console
uv sync
```

## Primary Development Commands

To check and resolve formatting issues in the codebase, run:

```console
uv run ruff format
```

To check and resolve linting issues in the codebase, run:

```console
uv run ruff check --fix
```

To check the unit tests in the codebase, run:

```console
uv run pytest
```

To check the typing in the codebase, run:

```console
uv run mypy
```

To generate a code coverage report after testing locally, run:

```console
uv run coverage html
```

To check the lock file is up to date:

```console
uv lock --check
```

## Shortcut Task Commands

###### For Running Individual Checks

```console
uv run poe check-lock
uv run poe check-format
uv run poe check-lint
uv run poe check-tests
uv run poe check-typing
```

###### For Running All Checks

```console
uv run poe check-all
```

###### For Running Individual Fixes

```console
uv run poe fix-format
uv run poe fix-lint
```

###### For Running All Fixes

```console
uv run poe fix-all
```

###### For Running All Fixes and Checks

```console
uv run poe fix-and-check-all
```

## Running the Exome Performance Benchmark

A benchmark is included which loads all known genes from hg38 into the `OverlapDetector`.
Each gene is then tested against the detector.
To run the benchmark:

```console
uv run pytest --benchmark
```

## Creating a release on PyPI

> [!NOTE]
> This project follows [Semantic Versioning](https://semver.org/), aka SemVer. In brief:
> 
> - MAJOR version when you make incompatible API changes
> - MINOR version when you add functionality in a backwards compatible manner
> - PATCH version when you make backwards compatible bug fixes

> [!IMPORTANT]
> Consider editing the changelog if there are any errors or necessary enhancements.

1. Clone the repository recursively, ensure you are on the main branch, and that the working directory is clean.
2. Check out a new branch to prepare the library for release.
3. Bump the version of the library to the desired SemVer using [`uv version`](https://docs.astral.sh/uv/reference/cli/#uv-version).
4. Commit the version bump changes with a Git commit message like `chore(release): bump to #.#.#`.
5. Push the commit, open a PR, ensure tests pass, and seek reviews.
6. Squash merge the PR into the `main` branch.
7. Tag the new commit on the main branch with the bumped version number.

> [!WARNING]
> The tag **must** be a valid SemVer version number and **must** match the version set by `uv run hatch version` in (3). The [publishing GitHub Action](.github/workflows/publish_pybedlite) is activated by a new tag on the `main` branch containing a valid SemVer version.

GitHub Actions will take care of the remainder of the deployment and release process:

1. Unit tests will be re-run.
2. A source distribution will be built.
3. A multi-arch multi-Python wheel (binary) distributions will be built.
4. Assets will be deployed to PyPi with the new version.
5. A [Conventional Commit](https://www.conventionalcommits.org/en/v1.0.0/)-aware changelog will be drafted.
6. A GitHub release will be created with the new version tag and the drafted changelog.

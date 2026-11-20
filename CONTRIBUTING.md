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
> The tag **must** be a valid SemVer version number (`X.Y.Z`, no `v` prefix), **must** match the version set by `uv version` in (3), and **must** point to a commit on `main`. The [publishing GitHub Action](.github/workflows/publish_pybedlite.yml) runs when such a tag is pushed, and fails (rather than skipping) when the tag is not on `main` or does not match the package version. It can also be started manually to retry a release (see below).

GitHub Actions will take care of the remainder of the deployment and release process:

1. The tag will be validated against `main` and the package version.
2. Unit tests will be re-run.
3. A source distribution will be built.
4. A multi-arch multi-Python wheel (binary) distributions will be built.
5. Assets will be deployed to PyPi with the new version.
6. A [Conventional Commit](https://www.conventionalcommits.org/en/v1.0.0/)-aware changelog will be drafted.
7. A GitHub release will be created with the new version tag and the drafted changelog.

### Retrying a release after a workflow fix

Merge the workflow fix into `main`, then run the `publish` workflow from `main` with the existing release tag:

```console
gh workflow run publish_pybedlite.yml --ref main -f tag=<X.Y.Z>
```

The workflow validates that the tag is on `main` and matches the package version, then tests and builds the tagged source using the workflow definitions (Python matrix, tool versions, docs build) from `main`. The tag and package version are not changed. Because `main`'s CI runs against the tagged code, retries are intended for recent releases whose code still passes `main`'s checks.

Existing PyPI files are skipped, so a partially published release can be retried. If a GitHub release for the tag already exists, its notes are replaced with a freshly generated changelog. Only the newest release tag is marked as the repository's latest GitHub release, so retrying an older tag does not demote a newer one.


[![Language][language-badge]][language-link]
[![Code Style][code-style-badge]][code-style-link]
[![Type Checked][type-checking-badge]][type-checking-link]
[![PEP8][pep-8-badge]][pep-8-link]
[![Code Coverage][code-coverage-badge]][code-coverage-link]
[![License][license-badge]][license-link]

---

[![Python package][python-package-badge]][python-package-link]
[![PyPI version][pypi-badge]][pypi-link]
[![PyPI download total][pypi-downloads-badge]][pypi-downloads-link]

---

[language-badge]:       http://img.shields.io/badge/language-python-brightgreen.svg
[language-link]:        http://www.python.org/
[code-style-badge]:     https://img.shields.io/badge/code%20style-black-000000.svg
[code-style-link]:      https://black.readthedocs.io/en/stable/ 
[type-checking-badge]:  http://www.mypy-lang.org/static/mypy_badge.svg
[type-checking-link]:   http://mypy-lang.org/
[pep-8-badge]:          https://img.shields.io/badge/code%20style-pep8-brightgreen.svg
[pep-8-link]:           https://www.python.org/dev/peps/pep-0008/
[code-coverage-badge]:  https://codecov.io/gh/fulcrumgenomics/pybedlite/branch/master/graph/badge.svg
[code-coverage-link]:   https://codecov.io/gh/fulcrumgenomics/pybedlite
[license-badge]:        http://img.shields.io/badge/license-MIT-blue.svg
[license-link]:         https://github.com/fulcrumgenomics/pybedlite/blob/master/LICENSE
[python-package-badge]: https://github.com/fulcrumgenomics/pybedlite/workflows/Python%20package/badge.svg
[python-package-link]:  https://github.com/fulcrumgenomics/pybedlite/actions?query=workflow%3A%22Python+package%22
[pypi-badge]:           https://badge.fury.io/py/pybedlite.svg
[pypi-link]:            https://pypi.python.org/pypi/pybedlite
[pypi-downloads-badge]: https://img.shields.io/pypi/dm/pybedlite
[pypi-downloads-link]:  https://pypi.python.org/pypi/pybedlite

# pybedlite

`pip install pybedlite`

**Requires python 3.6+**

# Getting Setup

[Poetry][poetry-link] is used to manage the python development environment. 

A simple way to create an environment with the desired version of python and poetry is to use [conda][conda-link].  E.g.:

```bash
conda create -n pybedlite python=3.6 poetry
conda activate pybedlite
poetry install
```

If, during `poetry install` on Mac OS X errors are encountered running gcc/clang to build `pybedtools` or other packages with native code, try setting the following and re-running `poetry install`:
```bash
export CFLAGS="-stdlib=libc++"
``` 

[poetry-link]: https://github.com/python-poetry/poetry
[conda-link]:  https://docs.conda.io/en/latest/miniconda.html

## Checking the Build
./ci/check.sh
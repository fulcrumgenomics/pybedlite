"""
Import shim for :mod:`superintervals`.

``superintervals`` ships a C++ extension built from source. Some interpreters
report ``CXX=gcc`` in ``sysconfig`` (notably Debian's ``python:3.14-slim``
images), which links the extension without ``libstdc++``. Importing it then
fails with ``undefined symbol: _ZTVN10__cxxabiv117__class_type_infoE``. Loading
``libstdc++`` globally before the import satisfies that symbol.

Import ``IntervalSet`` from here rather than from ``superintervals`` directly so
this preload always runs first.
"""

import ctypes
import platform

if platform.system() == "Linux":
    try:
        ctypes.CDLL("libstdc++.so.6", mode=ctypes.RTLD_GLOBAL)
    except OSError:
        pass

from superintervals import IntervalSet

__all__ = ["IntervalSet"]

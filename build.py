from setuptools import Extension
from Cython.Build import build_ext

ext_modules = [
    Extension(
        "cgranges",
        sources=["cgranges/python/cgranges.pyx", "cgranges/cgranges.c"],
        depends=[
            "cgranges/cgranges.h",
            "cgranges/khash.h",
            "cgranges/python/cgranges.pyx"
        ],
        include_dirs=["cgranges"]
    )
]

def build(setup_kwargs):
    """
    This function is mandatory in order to build the extensions.
    """
    setup_kwargs.update(
        {
            "ext_modules": ext_modules,
            "cmdclass": {"build_ext": build_ext},
        }
    )
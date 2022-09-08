from setuptools import setup, Extension
from Cython.Build import cythonize

setup(
        name="alpha_complex_periodic",
        packages="alpha_complex_periodic",
        ext_modules=cythonize(
            Extension(
                "alpha_complex_periodic",
                sources=["alpha_complex_periodic.pyx", "alpha_complex_periodic_persistence.cpp"],
                language="c++",
                libraries=["alpha_complex_periodic", "gmp", "mpfr"],
                extra_link_args=["-L./"]
            )
        )
)

from setuptools import setup, Extension
from Cython.Build import cythonize

setup(
        name="alpha_complex_periodic",
        version="0.1.0",
        packages=["alpha_complex_periodic"],
        author="Aaron Ouellette",
        url="https://github.com/ajouellette/alpha_complex_periodic",
        description="Python wrapper for the 3D alpha complex from GUDHI",

        setup_requires=["Cython>=0.29"],
        ext_modules=cythonize(
            Extension(
                "alpha_complex_periodic",
                sources=["alpha_complex_periodic/alpha_complex_periodic.pyx",
                    "alpha_complex_periodic/alpha_complex_periodic_persistence.cpp"],
                language="c++",
                libraries=["gmp", "mpfr"],
                extra_compile_args=["-O3", "-march=native"]
            ),
            annotate=True
        )
)

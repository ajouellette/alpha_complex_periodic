# cython: language_level=3
from libcpp.vector cimport vector
from libcpp cimport bool


cdef extern from "alpha_complex_periodic_persistence.h":
    cdef vector[vector[vector[double]]] _calc_persistence "calc_persistence"(vector[vector[double]], int, double, bool, bool, double) nogil except +


def calc_persistence(vector[vector[double]] coords, int coeff_field=2, double min_persistence=0, precision="safe", boxsize=None):
    cdef bool fast
    cdef bool exact
    cdef double _boxsize = 0.

    if boxsize is not None:
        _boxsize = boxsize

    if precision == "safe":
        fast = False
        exact = False
    elif precision == "fast":
        fast = True
        exact = False
    elif precision == "exact":
        fast = False
        exact = True
    else:
        raise ValueError("Invalid value for argument precision")

    return _calc_persistence(coords, coeff_field, min_persistence, fast, exact, _boxsize)

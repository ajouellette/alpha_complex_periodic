# cython: language_level=3
from libcpp.vector cimport vector


cdef extern from "alpha_complex_periodic_persistence.h":
    cdef vector[vector[vector[double]]] _calc_persistence "calc_persistence"(vector[vector[double]], int, double, double) except +


def calc_persistence(vector[vector[double]] coords, int coeff_field=2, double min_persistence=0, precision="safe", boxsize=None):
    if boxsize is None:
        boxsize = 0

    if precision == "safe":
        fast = False
        exact = False
    elif precision == "fast":
        fast = True
        exact = False
        raise NotImplementedError()
    elif precision == "exact":
        fast = False
        exact = True
        raise NotImplementedError()
    else:
        raise ValueError("Invalid value for argument precision")

    return _calc_persistence(coords, coeff_field, min_persistence, boxsize)

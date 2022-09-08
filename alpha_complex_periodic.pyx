# cython: language_level=3
from libcpp.vector cimport vector


cdef extern from "alpha_complex_periodic_persistence.h":
    cdef vector[vector[vector[double]]] _calc_persistence "calc_persistence"(vector[vector[double]], int, double, double) except +


def calc_persistence(vector[vector[double]] coords, int coeff_field=2, double min_persistence=0, boxsize=None):
    if boxsize is None:
        boxsize = 0

    return _calc_persistence(coords, coeff_field, min_persistence, boxsize)

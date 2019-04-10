"""
multivar.py

Multiplying multivariate polynomials
Kat Cannon-MacMartin | 4/10/19 | Marlboro College
"""

from numpy import fft

import time

def fexpand(pn):
    ret_pn = pn[0]
    for factor in xrange(1, len(pn)):
        ret_pn = fft.ifft(fft.fft(ret_pn) * fft.fft(pn[factor]))
    return ret_pn

def array_to_sage(co_array):
    sage_poly = 1
    for factor in co_array:
        sage_factor = 0
        for i in xrange(len(factor)):
            sage_factor += factor[i]*x**i
        sage_poly = sage_poly * sage_factor
    return sage_poly

def sage_to_array(sagepoly):
    sage_vars = list(sagepoly.variables())
    sage_fac_list = sagepoly[0].factor_list()
    

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
        

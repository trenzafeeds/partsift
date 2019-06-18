"""
threads.py
"""

from sage.symbolic.expression_conversions import polynomial as type_change
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import RationalField

import sys

from multiprocessing import Process, Queue, Pool

def strip(poly):
    return poly

def product(term1, term2, q):
    q.put(strip(term1 * term2))

def mult(term1, term2):
    q = Queue()
    p = Process(target=product, args=(term1, term2, q))
    p.start()
    print q.get()
    p.join()

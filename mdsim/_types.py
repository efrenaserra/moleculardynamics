# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 21:14:48 2019

@author: Efren A. Serra
"""
import math, numpy as np
from enum import Enum

__ALL__ = [
        'Mol' ,
        'Prop',
        'VecI',
        'VecR',
        'r_diff_vectorized',
        'rv_diff_vectorized',
        'ra_diff_vectorized',
        ]

class Prop(object):
    """
    Parameters
    ----------
        val  :
        sum  :
        sum2 : 
    """
    def __init__(self):
        self.val : float = 0.
        self.sum : float = 0.
        self.sum2: float = 0.
    
    """
    Parameters
    ----------
    """
    def accum(self):
        self.sum += self.val
        self.sum2 += math.sqrt(self.val)
    
    """
    Parameters
    ----------
    """
    def zero(self):
        self.sum = self.sum2 = 0.

    """
    Parameters
    ----------
        n    :
    """
    def avg(self, n: int):
        self.sum /= n
        self.sum2 = math.sqrt(max([self.sum2 / n - math.sqrt(self.sum), 0.]))
    
    """
    Parameters:
    """
    def est(self):
        return "%7.4f %7.4f"%(self.sum, self.sum2)

    """
    Parameters:
    """
    def __repr__(self):
        return "<val: %f; sum: %f; sum2: %f>"%(self.val,self.sum,self.sum2)

class VecI(object):
    def __init__(self, x : int=0, y : int=0):
        self.x = x
        self.y = y

    def __repr__(self):
        return "<x: %d, y: %d>"%(self.x,self.y)

class VecR(object):
    def __init__(self, x : float=0., y : float=0.):
        self.x = x
        self.y = y

    def __sub__(self, other):
        return VecR((self.x - other.x), (self.y - other.y))

    def vcsum(self):
        return self.x + self.y

    def __repr__(self):
        return "<x: %f, y: %f>"%(self.x,self.y)

class Mol(object):
    def __init__(self, r: VecR, rv: VecR, ra: VecR):
        self.r  = r
        self.rv = rv
        self.ra = ra

    def __repr__(self):
        return "<r: %s, rv: %s, ra: %s>"%(self.r,self.rv,self.ra)

def r_diff(a,b):
    """Return molecular position difference.
    """
    return VecR(a.r.x - b.r.x, a.r.y - b.r.y)

def rv_diff(a,b):
    """Return molecular velocity difference.
    """
    return VecR(a.rv.x - b.rv.x, a.rv.y - b.rv.y)

def ra_diff(a,b):
    """Return molecular acceleration difference.
    """
    return VecR(a.ra.x - b.ra.x, a.ra.y - b.ra.y)

r_diff_vectorized = np.vectorize(r_diff,otypes=[VecR],\
                                 signature='(),()->()')

rv_diff_vectorized = np.vectorize(rv_diff, otypes=[VecR],\
                                 signature='(),()->()')

ra_diff_vectorized = np.vectorize(ra_diff, otypes=[VecR],\
                                 signature='(),()->()')

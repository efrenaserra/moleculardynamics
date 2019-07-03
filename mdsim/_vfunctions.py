# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 10:53:02 2019

@author: Efren A. Serra
"""

import math, numpy as np, random
import _types
from _types import VecR, Mol

__ALL__ = [
        'r_diff_vectorized',
        'rv_add_vectorized',
        'rv_diff_vectorized',
        'rv_dot_vectorized',
        'rv_rand_vectorized',
        'rv_scale_vectorized',
        'rv_sadd_vectorized',
        'ra_diff_vectorized',
        'leapfrog_integrate_rv_vectorized'
        'leapfrog_integrate_r_vectorized'
        ]

def r_diff(a, b):
    """Return molecular position difference.
    """
    return VecR(a.r.x - b.r.x, a.r.y - b.r.y)

def r_wrap(m,region):
    """
    Parameters
    ----------
    m : Mol, 
    region : VecR, 
    """
    # wrap the x-coordinate
    if m.r.x >= 0.5 * region.x:
        m.r.x -= region.x
    elif m.r.x < -0.5 * region.x:
        m.r.x += region.x

    # wrap the y-coordinate
    if m.r.y >= 0.5 * region.y:
        m.r.y -= region.y
    elif m.r.y < -0.5 * region.y:
        m.r.y += region.y

def vecr_wrap(vecr,region):
    """
    Parameters
    ----------
    m : Mol, 
    region : VecR, 
    """
    # wrap the x-coordinate
    if vecr.x >= 0.5 * region.x:
        vecr.x -= region.x
    elif vecr.x < -0.5 * region.x:
        vecr.x += region.x

    # wrap the y-coordinate
    if vecr.y >= 0.5 * region.y:
        vecr.y -= region.y
    elif vecr.y < -0.5 * region.y:
        vecr.y += region.y

    return vecr

def rv_diff(a,b):
    """Return molecular velocity difference.
    """
    return VecR(a.rv.x - b.rv.x, a.rv.y - b.rv.y)

def ra_diff(a,b):
    """Return molecular acceleration difference.
    """
    return VecR(a.ra.x - b.ra.x, a.ra.y - b.ra.y)

def rv_rand(m):
    """Set molecular velocity components.
    """
    s : float = 2. * math.pi * random.random()
    m.rv.x = math.cos(s)
    m.rv.y = math.sin(s)

def rv_scale(m, s):
    """Set molecular velocity components.
    """
    m.rv.x *= s
    m.rv.y *= s

def rv_add(v, m):
    """Accumulate molecular velocity components.
    """
    v.x += m.rv.x
    v.y += m.rv.y

def rv_dot(a, b):
    """Accumulate molecular velocity components.
    """
    v : float = 0.
    for i in range(a.shape[0]):
        v += (a[i].rv.x * b[i].rv.x + a[i].rv.y * b[i].rv.y)
    return v

def rv_sadd(array, s, v):
    """Scale molecular velocity components.
    """
    for m in array:
        m.rv.x += (s * v.x)
        m.rv.y += (s * v.y)

def ra_zero(m):
    """Zero molecular acceleration components.
    """
    m.ra.x = m.ra.y = 0.

def leapfrog_integrate_rv(m, s):
        m.rv.x += (s * m.ra.x)
        m.rv.y += (s * m.ra.y)
leapfrog_integrate_rv_vectorized = np.vectorize(leapfrog_integrate_rv)

def leapfrog_integrate_r(m, s):
        m.r.x += (s * m.rv.x)
        m.r.y += (s * m.rv.y)
leapfrog_integrate_r_vectorized = np.vectorize(leapfrog_integrate_r)

r_diff_vectorized = np.vectorize(r_diff)
r_wrap_vectorized = np.vectorize(r_wrap)

rv_add_vectorized  = np.vectorize(rv_add)
rv_diff_vectorized = np.vectorize(rv_diff, otypes=[VecR],\
                                 signature='(),()->()')
rv_dot_vectorized   = np.vectorize(rv_dot, otypes=[float], signature='(n),(n)->()')
rv_rand_vectorized  = np.vectorize(rv_rand, signature='()->()')
rv_scale_vectorized = np.vectorize(rv_scale, signature='(),()->()')
rv_sadd_vectorized  = np.vectorize(rv_sadd, signature='(n),(),()->()')

ra_diff_vectorized = np.vectorize(ra_diff, otypes=[VecR],\
                                 signature='(),()->()')
ra_zero_vectorized = np.vectorize(ra_zero)
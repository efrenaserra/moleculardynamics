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
        'leapfrog_integrate_rv_vectorized',
        'leapfrog_integrate_r_vectorized',
        'vecr_div',
        'vecr_dot',
        'vecr_mul',
        'vecr_sadd',
        'vecr_wrap',
        ]

def r_diff(a, b):
    """Return molecular position difference.
    """
    return VecR(a.r.x - b.r.x, a.r.y - b.r.y)

r_diff_vectorized = np.vectorize(r_diff)

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

r_wrap_vectorized = np.vectorize(r_wrap)

def vecr_sadd(a, s, v):
    """Scale molecular velocity components.
    """
    a.x += (s * v.x)
    a.y += (s * v.y)

    return a

def vecr_dot(a,b):
    """Divide two VecR objects component-wise.
    Parameters
    ----------
    a : VecR, a molecule position
    b : VecR, another molecule position
    """
    return (a.x * b.x + a.y * b.y)

def vecr_div(a,b):
    """Divide two VecR objects component-wise.
    Parameters
    ----------
    a : VecR, a molecule position
    b : VecR, another molecule position
    """
    return VecR(a.x / b.x, a.y / b.y)

def vecr_mul(a,b):
    """Multiply two VecR objects component-wise.
    Parameters
    ----------
    a : VecR, a molecule position
    b : VecR, another molecule position
    """
    return VecR(a.x * b.x, a.y * b.y)

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

rv_diff_vectorized = np.vectorize(rv_diff)

def ra_diff(a,b):
    """Return molecular acceleration difference.
    """
    return VecR(a.ra.x - b.ra.x, a.ra.y - b.ra.y)

def rv_rand(m):
    """Set molecular velocity components.
    Parameters
    ----------
    m : Mol, the molecular object
    """
    s : float = 2. * math.pi * random.random()
    m.rv = VecR(math.cos(s), math.sin(s))

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

rv_add_vectorized  = np.vectorize(rv_add)

def rv_dot(a, b):
    """Accumulate molecular velocity components.
    """
    return (a.rv.x * b.rv.x + a.rv.y * b.rv.y)

rv_dot_vectorized   = np.vectorize(rv_dot)

def rv_sadd(m, s, v):
    """Scale molecular velocity components.
    """
    m.rv.x += (s * v.x)
    m.rv.y += (s * v.y)

rv_sadd_vectorized  = np.vectorize(rv_sadd)

def ra_zero(m):
    """Zero molecular acceleration components.
    """
    m.ra = VecR()

def leapfrog_integrate_rv(m, s):
        m.rv.x += (s * m.ra.x)
        m.rv.y += (s * m.ra.y)
leapfrog_integrate_rv_vectorized = np.vectorize(leapfrog_integrate_rv)

def leapfrog_integrate_r(m, s):
        m.r.x += (s * m.rv.x)
        m.r.y += (s * m.rv.y)
leapfrog_integrate_r_vectorized = np.vectorize(leapfrog_integrate_r)


rv_rand_vectorized  = np.vectorize(rv_rand)
rv_scale_vectorized = np.vectorize(rv_scale)

ra_diff_vectorized = np.vectorize(ra_diff, otypes=[VecR],\
                                 signature='(),()->()')
ra_zero_vectorized = np.vectorize(ra_zero)
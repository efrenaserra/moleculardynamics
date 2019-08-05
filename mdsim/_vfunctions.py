# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 10:53:02 2019

@author: Efren A. Serra
"""

import math, numpy as np, random
from _types import VecR, Mol

__ALL__ = [
        'r_diff',
        'ra_sadd',
        'rv_add',
        'rv_rand',
        'rv_sadd',
        'rv_scale',
        'vecr_div',
        'vecr_dot',
        'vecr_mul',
        'vecr_sadd',
        'vecr_wrap',
        ]

def r_diff(a: Mol, b: Mol):
    """Return molecular relative vector difference.
    """
    return VecR(a.r.x - b.r.x, a.r.y - b.r.y)

def r_wrap(m: Mol, region: VecR):
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

def vecr_sadd(a: VecR, s: float, v: VecR):
    """Scale molecular velocity components.
    """
    a.x += (s * v.x)
    a.y += (s * v.y)

    return a

def vecr_dot(a: VecR, b: VecR):
    """Divide two VecR objects component-wise.
    Parameters
    ----------
    a : VecR, a molecule position
    b : VecR, another molecule position
    """
    return (a.x * b.x + a.y * b.y)

def vecr_div(a: VecR, b: VecR):
    """Divide two VecR objects component-wise.
    Parameters
    ----------
    a : VecR, a molecule position
    b : VecR, another molecule position
    """
    return VecR(a.x / b.x, a.y / b.y)

def vecr_mul(a: VecR, b: VecR):
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
    return (a.rv.x * b.rv.x + a.rv.y * b.rv.y)

def rv_sadd(m, s, v):
    """Scale molecular velocity components.
    """
    m.rv.x += (s * v.x)
    m.rv.y += (s * v.y)

def ra_sadd(m, s, v):
    """Scale molecular acceleration components.
    """
    m.ra.x += (s * v.x)
    m.ra.y += (s * v.y)

def ra_zero(m):
    """Zero molecular acceleration components.
    """
    m.ra.x = 0.0
    m.ra.y = 0.0

# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 10:53:02 2019

/*********************************************************************

  This program is copyright material accompanying the book
  "The Art of Molecular Dynamics Simulation", 2nd edition,
  by D. C. Rapaport, published by Cambridge University Press (2004).

  Copyright (C) 2004, 2011  D. C. Rapaport

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************/

@author: Efren A. Serra
"""

import math, random
from _types import VecR, Mol

__ALL__ = [
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
        'vrand',
        ]

def vrand(p):
    """Produce unit vectors in two dimensions.
    """
    s : float = 2. * math.pi * random.random()
    p.rv.x = math.cos(s)
    p.rv.y = math.sin(s)

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
    """Set molecular velocity components to random values.
    Parameters
    ----------
    m : Mol, the molecular object
    """
    s : float = 2. * math.pi * random.random()
    m.rv.x = math.cos(s)
    m.rv.y = math.sin(s)

def rv_scale(m, s):
    """Scale molecular velocity components.
    """
    m.rv.x *= s
    m.rv.y *= s

def rv_add(v, m):
    """Accumulate molecular velocity components.
    """
    v.x += m.rv.x
    v.y += m.rv.y

def rv_dot(a, b):
    """Lenght squared of velocity vector.
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

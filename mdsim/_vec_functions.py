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

from _types import VecI, VecR, Mol

__ALL__ = [
        'ra_sadd',
        'rv_add',
        'rv_sadd',
        'rv_scale',
        'vecr_dot',
        'vecr_sadd',
        'vecr_wrap',
        ]

def r_wrap(m: Mol, region: VecR):
    """
    Parameters
    ----------
    m : Mol, 
    region : VecR, 
    """
    # Wrap the x-coordinate
    if m.r.x >= 0.5 * region.x:
        m.r.x -= region.x
    elif m.r.x < -0.5 * region.x:
        m.r.x += region.x

    # Wrap the y-coordinate
    if m.r.y >= 0.5 * region.y:
        m.r.y -= region.y
    elif m.r.y < -0.5 * region.y:
        m.r.y += region.y

    # Wrap the z-coordinate
    if m.r.z >= 0.5 * region.z:
        m.r.z -= region.z
    elif m.r.z < -0.5 * region.z:
        m.r.z += region.z

def _VCell_wrap_all(vc: VecI, cells: VecI, rs: VecR, region: VecR):
    # Wrap the x-coordinate
    if vc.x >= cells.x:
        vc.x = 0
        rs.x = region.x
    elif vc.x < 0:
        vc.x = cells.x - 1
        rs.x = - region.x

    # Wrap the y-coordinate
    if vc.y >= cells.y:
        vc.y = 0
        rs.y = region.y
    elif vc.y < 0:
        vc.y = cells.y - 1
        rs.y = - region.y

    # Wrap the z-coordinate
    if vc.z >= cells.z:
        vc.z = 0
        rs.z = region.z
    elif vc.z < 0:
        vc.z = cells.z - 1
        rs.z = - region.z

def veci_mul(a: VecR, b: VecR):
    """Multiply two VecR objects component-wise.
    Parameters
    ----------
    a : VecR, a molecule position
    b : VecR, another molecule position
    """
    return VecI(int(a.x * b.x), int(a.y * b.y), int(a.z * b.z))

def vecr_sadd(a: VecR, s: float, v: VecR):
    """Scale molecular velocity components.
    """
    a.x += (s * v.x)
    a.y += (s * v.y)
    a.z += (s * v.z)

    return a

def vecr_dot(a: VecR, b: VecR):
    """Divide two VecR objects component-wise.
    Parameters
    ----------
    a : VecR, a molecule position
    b : VecR, another molecule position
    """
    return (a.x * b.x + a.y * b.y + a.z * b.z)

def vecr_wrap(vecr,region):
    """
    Parameters
    ----------
    m : Mol, 
    region : VecR, 
    """
    # Wrap the x-coordinate
    if vecr.x >= 0.5 * region.x:
        vecr.x -= region.x
    elif vecr.x < -0.5 * region.x:
        vecr.x += region.x

    # Wrap the y-coordinate
    if vecr.y >= 0.5 * region.y:
        vecr.y -= region.y
    elif vecr.y < -0.5 * region.y:
        vecr.y += region.y

    # Wrap the z-coordinate
    if vecr.z >= 0.5 * region.z:
        vecr.z -= region.z
    elif vecr.z < -0.5 * region.z:
        vecr.z += region.z

    return vecr

def rv_diff(a,b):
    """Return molecular velocity difference.
    """
    return VecR(a.rv.x - b.rv.x, a.rv.y - b.rv.y, a.rv.z - b.rv.z)

def ra_diff(a,b):
    """Return molecular acceleration difference.
    """
    return VecR(a.ra.x - b.ra.x, a.ra.y - b.ra.y, a.ra.z - b.ra.z)

def rv_scale(m, s):
    """Scale molecular velocity components.
    """
    m.rv.x *= s
    m.rv.y *= s
    m.rv.z *= s

def rv_add(v, m):
    """Accumulate molecular velocity components.
    """
    v.x += m.rv.x
    v.y += m.rv.y
    v.z += m.rv.z

def rv_dot(a, b):
    """Lenght squared of velocity vector.
    """
    return (a.rv.x * b.rv.x + a.rv.y * b.rv.y + a.rv.z * b.rv.z)

def rv_sadd(m, s, v):
    """Scale molecular velocity components.
    """
    m.rv.x += (s * v.x)
    m.rv.y += (s * v.y)
    m.rv.z += (s * v.z)

def ra_sadd(m, s, v):
    """Scale molecular acceleration components.
    """
    m.ra.x += (s * v.x)
    m.ra.y += (s * v.y)
    m.ra.z += (s * v.z)

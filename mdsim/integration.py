# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 11:04:01 2019

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
from .core.types import RVec

__all__ = [
        'leapfrog_update_coordindates',
        'leapfrog_update_velocities',
        ]

def leapfrog_update_coordinates(r: RVec, scale: float, rv: RVec):
    """Integrate the coordinates using the Leapfrog scheme."""
    r.x += (scale * rv.x)
    r.y += (scale * rv.y)
    r.z += (scale * rv.z)

def leapfrog_update_velocities(rv: RVec, scale: float, ra: RVec):
    """Integrate the coordinates using the Leapfrog scheme."""
    rv.x += (scale * ra.x)
    rv.y += (scale * ra.y)
    rv.z += (scale * ra.z)
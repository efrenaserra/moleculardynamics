# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:51:08 2019

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

@author: serra
"""

import math, random

__ALL__ = [
        'rv_rand',
        ]


def rv_rand(m):
    """Set molecular velocity components to random values.
    Parameters
    ----------
    m : Mol, the molecular object
    """
    s : float = 2.
    x : float = 0.
    y : float = 0.
    while s > 1.:
        x = 2. * random.random() - 1
        y = 2. * random.random() - 1
        s = x * x + y * y

    m.rv.z = 1. - 2. * s
    s = 2. * math.sqrt(1. - s)
    m.rv.x = s * x
    m.rv.y = s * y

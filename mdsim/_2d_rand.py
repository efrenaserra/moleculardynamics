# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:57:55 2019

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

@author: Efren Antonio Serra (Read)
"""

import math, random
from .core.types import Molecule

__all__ = ['rv_rand']


def rv_rand(m: Molecule) -> None:
    """Set molecular velocity components to random values.
    Parameters
    ----------
    m : Mol, the molecular object
    """
    s : float = 2. * math.pi * random.random()
    m.rv.x = math.cos(s)
    m.rv.y = math.sin(s)
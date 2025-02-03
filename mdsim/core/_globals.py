# -*- coding: utf-8 -*-
"""
Module defining global variables.

Created on Tue Jun 25 14:31:17 2019

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
from mdsim.core.types import VecI
# Disallow reloading this module so as to prevent re-initializing these global
# variables
if '_is_loaded' in globals():
    raise RuntimeError('Reloading mdsim._globals is not allowed')
_is_loaded = True


"""Global variables
"""
_mdsim_globals = {}
_namelist_converter = {
    'deltaT'         : lambda x: float(x),
    'density'        : lambda x: float(x),
    'initUcell'      : lambda x,y,z: VecI(int(x),int(y),int(z)),
    'limitVel'       : lambda x: int(x),
    'nebrTabFac'     : lambda x: int(x),
    'randSeed'       : lambda x: int(x),
    'rangeVel'       : lambda x: float(x),
    'rNebrShell'     : lambda x: float(x),
    'sizeHistVel'    : lambda x: int(x),
    'stepAvg'        : lambda x: int(x),
    'stepEquil'      : lambda x: int(x),
    'stepInitlzTemp' : lambda x: int(x),
    'stepLimit'      : lambda x: int(x),
    'stepVel'        : lambda x: int(x),
    'temperature'    : lambda x: float(x),
    }
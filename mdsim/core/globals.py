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
from __future__ import annotations
from typing import Callable, TypedDict

from .types import IVec
# Disallow reloading this module so as to prevent re-initializing these global
# variables
if '_is_loaded' in globals():
    raise RuntimeError('Reloading mdsim._globals is not allowed')
_is_loaded = True


"""Global variables
"""
#_mdsim_globals = {}
class ConverterDict(TypedDict):
    delta_t: Callable[[float], float]
    density: Callable[[float], float]
    init_U_cell: Callable[[float, float, float], IVec]
    limit_vel: Callable[[float], int]
    nbr_tabfac: Callable[[float], int]
    rand_seed: Callable[[float], int]
    range_vel: Callable[[float], float]
    rnbr_shell: Callable[[float], float]
    sizeHist_vel: Callable[[float], int]
    step_avg: Callable[[float], int]
    step_equil: Callable[[float], int]
    setp_itemp: Callable[[float], int]
    step_limit: Callable[[float], int]
    step_vel: Callable[[float], int]
    temperature: Callable[[float], float]

_namelist_converter: ConverterDict = {
    'delta_t'        : lambda x: float(x),
    'density'        : lambda x: float(x),
    'init_U_cell'    : lambda x,y,z: IVec(x,y,z),
    'limit_vel'      : lambda x: int(x),
    'nbr_TabFac'     : lambda x: int(x),
    'rand_seed'      : lambda x: int(x),
    'range_vel'      : lambda x: float(x),
    'rNebr_shell'    : lambda x: float(x),
    'sizeHist_vel'   : lambda x: int(x),
    'step_avg'       : lambda x: int(x),
    'step_equil'     : lambda x: int(x),
    'stepInitlzTemp' : lambda x: int(x),
    'step_limit'     : lambda x: int(x),
    'step_vel'       : lambda x: int(x),
    'temperature'    : lambda x: float(x),
    }
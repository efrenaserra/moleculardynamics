"""A molecular dynamics simulator module.
"""

__all__ = [
    'physics',
    ]

import logging
from typing import Dict

logger = logging.getLogger("mdsim.__init__")
from . import physics
from .core.types import IVec

# Dictionary/Registry of Global Variables
_MDSIM_GLOBALS: Dict = {}
# Dictionary of Conversters
_NAMELIST_CONVERTER: Dict = {
    'delta_t'        : lambda x: float(x),
    'density'        : lambda x: float(x),
    'init_u_cell'    : lambda x,y,z: IVec(int(x),int(y),int(z)),
    'limit_vel'      : lambda x: int(x),
    'nebrTabFac'     : lambda x: int(x),
    'randSeed'       : lambda x: int(x),
    'rangeVel'       : lambda x: float(x),
    'rNebrShell'     : lambda x: float(x),
    'sizeHistVel'    : lambda x: int(x),
    'step_avg'       : lambda x: int(x),
    'step_equil'     : lambda x: int(x),
    'stepInitlzTemp' : lambda x: int(x),
    'step_limit'     : lambda x: int(x),
    'step_vel'       : lambda x: int(x),
    'temperature'    : lambda x: float(x),
    }
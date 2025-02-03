from __future__ import annotations
from typing import TypedDict

import re

class MDUniverse(object):
    """
    The MD simulation simple system class.
    """
    __slots__ = ()

    def __init__(self):
        """
        """
        self._variable_registry = self.MDRegistry({})

    class MDRegistry(TypedDict):
        delta_t: float
        density: float
        init_U_cell: tuple[int, int] # The unit cell
        step_avg: int
        step_equil: int
        step_limit: int
        temperature: float

    def __getattr__(self, name):
        """
        """
        pass

    def __repr__(self):
        """
        """
        pass

    def setup(self, inparams):
        """
        Parameters
        ==========
        inparams : str, the filename
        """
        with open(inparams, 'r') as f:
            pattern = re.compile(r'init_U_cell')
            for line in f:
                m = pattern.match(line)
                if not m:
                    (k, v) = line.split()
                    self._variable_registry[k] = _NAMELIST_CONVERTER[k](v)
                else:
                    k = 'init_U_cell'
                    line = line[len(k):]
                    (nx, ny) = line.split()
                    # Matrix of molecular unit cells
                    self._variable_registry[k] = _NAMELIST_CONVERTER[k](nx, ny)

    def single_step(self):
        """
        """
        pass

    def exercise(self):
        """
        """
        self.setup()
        self.single_step()
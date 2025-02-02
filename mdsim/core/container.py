from __future__ import decorations

import re
from typing import Dict

class MDContainer(object):
    """
    The MD simulation container class.
    """
    __slots__ = ()

    def __init__(self):
        """
        """
        self._variable_registry: Dict = Dict()

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
                    (k,v) = line.split()
                    self._variable_registry[k] = _NAMELIST_CONVERTER[k](v)
                else:
                    k = 'init_U_cell'
                    line = line[len(k):]
                    (nx, ny) = line.split()
                    # Matrix of molecular unit cells
                    self._mdsim_var_registry[k] = _NAMELIST_CONVERTER[k](nx, ny)

    def single_step(self):
        """
        """
        pass

    def exercise(self):
        """
        """
        self.setup()
        self.single_step()
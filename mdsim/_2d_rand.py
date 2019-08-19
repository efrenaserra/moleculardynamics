# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:57:55 2019

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
    s : float = 2. * math.pi * random.random()
    m.rv.x = math.cos(s)
    m.rv.y = math.sin(s)

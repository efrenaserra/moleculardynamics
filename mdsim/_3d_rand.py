# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:51:08 2019

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

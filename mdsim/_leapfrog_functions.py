# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 11:04:01 2019

@author: Efren A. Serra
"""
import math, random
from _types import VecR

__ALL__ = [
        'vrand',
        'leapfrog_update_coordindates',
        'leapfrog_update_velocities',
        ]

def vrand(p):
    """Produce unit vectors in two dimensions.
    """
    s : float = 2. * math.pi * random.random()
    p.rv.x = math.cos(s)
    p.rv.y = math.sin(s)

def leapfrog_update_coordinates(r: VecR, scale: float, rv: VecR):
    """Integrate the coordinates using the Leapfrog scheme.
    """
    r.x += (scale * rv.x)
    r.y += (scale * rv.y)

def leapfrog_update_velocities(rv: VecR, scale: float, ra: VecR):
    """Integrate the coordinates using the Leapfrog scheme.
    """
    rv.x += (scale * ra.x)
    rv.y += (scale * ra.y)
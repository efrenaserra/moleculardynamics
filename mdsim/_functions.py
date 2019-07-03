# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 11:04:01 2019

@author: Efren A. Serra
"""
import math, random
from _types import VecR

def vrand(p):
    """Produce unit vectors in two dimensions.
    """
    s : float = 2. * math.pi * random.random()
    p.rv.x = math.cos(s)
    p.rv.y = math.sin(s)
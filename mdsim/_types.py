# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 21:14:48 2019

@author: Efren A. Serra
"""
import math

__ALL__ = [
        'Mol' ,
        'Prop',
        'VecI',
        'VecR',
        ]

class Prop(object):
    def __init__(self):
        """
        Parameters
        ----------
            val  :
            sum  :
            sum2 : 
        """
        self.val : float = 0.
        self.sum : float = 0.
        self.sum2: float = 0.
    
    def accum(self):
        """
        Parameters
        ----------
        """
        self.sum += self.val
        self.sum2 += math.sqrt(self.val)
    
    def zero(self):
        """
        Parameters
        ----------
        """
        self.sum = self.sum2 = 0.

    def avg(self, n: int):
        """
        Parameters
        ----------
            n    :
        """
        self.sum /= n
        self.sum2 = math.sqrt(max([self.sum2 / n - math.sqrt(self.sum), 0.]))
    
    def est(self):
        """
        Parameters:
        """
        return "%7.4f %7.4f"%(self.sum, self.sum2)

    def __repr__(self):
        """
        Parameters:
        """
        return "<val: %f; sum: %f; sum2: %f>"%(self.val,self.sum,self.sum2)

class VecI(object):
    def __init__(self, x : int=0, y : int=0):
        self.x = x
        self.y = y

    def __repr__(self):
        return "<x: %d, y: %d>"%(self.x,self.y)

class VecR(object):
    def __init__(self, x : float=0., y : float=0.):
        self.x = x
        self.y = y

    def __sub__(self, other):
        return VecR((self.x - other.x), (self.y - other.y))

    def vcsum(self):
        return self.x + self.y

    def zero(self):
        """
        Parameters
        ----------
        """
        self.x = self.y = 0.

    def __repr__(self):
        return "<x: %f, y: %f>"%(self.x,self.y)

class Mol(object):
    def __init__(self, r: VecR, rv: VecR, ra: VecR):
        self.r  = r
        self.rv = rv
        self.ra = ra

    def __repr__(self):
        return "<r: %s, rv: %s, ra: %s>"%(self.r,self.rv,self.ra)

    def __sub__(self, other):
        return VecR((self.r.x - other.r.x), (self.r.y - other.r.y))

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
            n:
        """
        self.sum /= n
        self.sum2 = math.sqrt(max([self.sum2 /  n - (self.sum * self.sum), 0.]))

    def est(self):
        """
        Parameters:
        """
        return self.sum, self.sum2

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

    def __sub__(self, other):
        """Rerturn relative vector difference
        """
        return VecI((self.x - other.x), (self.y - other.y))

class VecR(object):
    def __init__(self, x : float=0., y : float=0.):
        self.x = x
        self.y = y

    def __sub__(self, other):
        """Rerturn relative vector difference.
        """
        return VecR((self.x - other.x), (self.y - other.y))

    def vcsum(self):
        return self.x + self.y

    def zero(self):
        """Zero the vector components.
        """
        self.x = 0.
        self.y = 0.

    def __repr__(self):
        return "<x: %f, y: %f>"%(self.x,self.y)

class Mol(object):
    def __init__(self, r: VecR=None, rv: VecR=None, ra: VecR=None):
        self.r  = r
        self.rv = rv
        self.ra = ra

    def __repr__(self):
        return "<r: %s, rv: %s, ra: %s>"%(self.r,self.rv,self.ra)

    def r_diff(self, other):
        """Rerturn molecule's relative vector difference.
        """
        return self.r - other.r

    def r_wrap(self, region):
        """Rerturn molecule's relative vector difference.
        Parameters
        ----------
        region : VecR, 
        """
        # wrap the x-coordinate
        if self.r.x >= 0.5 * region.x:
            self.r.x -= region.x
        elif self.r.x < -0.5 * region.x:
            self.r.x += region.x

        # wrap the y-coordinate
        if self.r.y >= 0.5 * region.y:
            self.r.y -= region.y
        elif self.r.y < -0.5 * region.y:
            self.r.y += region.y

    def ra_zero(self):
        self.ra.zero()

    def update_coordinates(self, integration_scheme, *args):
        """Integrate the coordinates using scheme.
        """
        integration_scheme(self.r, *args, self.rv)

    def update_velocities(self, integration_scheme, *args):
        """Integrate the velocities using Leapfrog scheme.
        """
        integration_scheme(self.rv, *args, self.ra)
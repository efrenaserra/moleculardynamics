# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 21:14:48 2019

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
import math

__ALL__ = [
        'Mol' ,
        'Prop',
        'VecI',
        'VecR',
        ]

class Prop(object):
    """
    A class representing a prognostic variable, such as kinetic Energy,
    total Enery and temperature, in an MD simulation.

    Attributes
    ----------
    val  : float
        the property value
    sum  : float
        the sum of the property
    sum2 : float
        the sum squared of the cummulative property

    Methods
    -------
    accum()
        Accumulates value and value squared of property.
    zero()
        Zeros the sum and sum squared values of property.
    avg(n)
        Computes the average and variance of the property.
    """
    def __init__(self, val: float = 0.0, sum: float = 0.0, sum2: float = 0.0):
        """
        Parameters
        ----------
        val  : float
            the property value
        sum  : float
            the sum of the property
        sum2 : float
            the sum squared of the cummulative property
        """
        self.val : float = val
        self.sum : float = sum
        self.sum2: float = sum2
    
    def accum(self):
        """Accumulates value and value squared of property.
        """
        self.sum += self.val
        self.sum2 += self.val * self.val

    def zero(self):
        """Zeros the sum and sum squared values of property.
        """
        self.sum = self.sum2 = 0.

    def avg(self, n: int):
        """Computes the average and variance of the property.
        """
        self.sum /= n
        self.sum2 = math.sqrt(max([self.sum2 /  n - (self.sum * self.sum), 0.]))

    def est(self):
        """Returns a tuple representing the sum and sum squared of the property.
        """
        return self.sum, self.sum2

    def __repr__(self):
        """Returns a formatted string representation (val, sum, sum2) of the property.
        """
        return "<val: %f; sum: %f; sum2: %f>"%(self.val,self.sum,self.sum2)

class VecI(object):
    """
    A class representing a physical property, such as kinetic Energy,
    total Enery and temperature, in an MD simulation.

    Attributes
    ----------
    x  : float
        the property value
    y  : float
        the sum of the property
    """
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
    """
    A class representing a point in a 2-Dimensional Cartessian coordinate
    system of reference during an MD simulation.

    Attributes
    ----------
    x  : float
        the x coordinate
    y  : float
        the y coordinate
    """
    def __init__(self, x : float=0., y : float=0.):
        self.x = x
        self.y = y

    def __sub__(self, other):
        """Return relative vector difference.
        """
        return VecR((self.x - other.x), (self.y - other.y))

    def vcsum(self):
        return self.x + self.y

    def wrap(self, region):
        """
        Parameters
        ----------
        region : VecR, 
        """
        # Wrap the x-coordinate
        if self.x >= 0.5 * region.x:
            self.x -= region.x
        elif self.x < -0.5 * region.x:
            self.x += region.x

        # Wrap the y-coordinate
        if self.y >= 0.5 * region.y:
            self.y -= region.y
        elif self.y < -0.5 * region.y:
            self.y += region.y

        return self

    def zero(self):
        """Zero the vector components.
        """
        self.x = 0.
        self.y = 0.

        return self

    def __repr__(self):
        return "<x: %f, y: %f>"%(self.x,self.y)

class Mol(object):
    """
    A class representing a molecule, i.e., position: r=(x,y), velocity: rv=(vx,vy),
    and acceleration: ra=(ax,ay), in an MD simulation.

    Attributes
    ----------
    r  : float
        the property value
    rv  : float
        the sum of the property
    ra : float
        the sum squared of the cummulative property

    Methods
    -------
    r_diff(other)
        Accumulates value and value squared of property.
    r_wrap(region)
        Zeros the sum and sum squared values of property.
    ra_zero()
        Computes the average and variance of the property.
    update_coordinates(integration_scheme, *args)
        Integrate the coordinates using integration scheme function.
    update_velocities(integration_scheme, *args)
        Integrate the velocities using integration scheme function.
    """
    def __init__(self, r: VecR=None, rv: VecR=None, ra: VecR=None):
        self.r  = r
        self.rv = rv
        self.ra = ra

    def __repr__(self):
        return "<r: %s, rv: %s, ra: %s>"%(self.r,self.rv,self.ra)

    def r_diff(self, other):
        """Return molecule's relative vector difference.
        """
        return self.r - other.r

    def r_wrap(self, region):
        """Return molecule's relative vector difference.
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

        return self

    def ra_zero(self):
        self.ra.zero()

        return self

    def ra_sadd(self, s, v):
        """Scale molecular acceleration components.
        """
        self.ra.x += (s * v.x)
        self.ra.y += (s * v.y)

        return self

    def update_coordinates(self, integration_scheme, *args):
        """Integrate the coordinates using scheme.
        """
        integration_scheme(self.r, *args, self.rv)

        return self

    def update_velocities(self, integration_scheme, *args):
        """Integrate the velocities using scheme.
        """
        integration_scheme(self.rv, *args, self.ra)

        return self
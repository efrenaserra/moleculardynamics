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

@author: Efren Antonio Serra
"""
from __future__ import annotations

from typing import Generic, TypeAlias, TypeVar
import math

T = TypeVar('T')
class Point(Generic[T]):
    def __init__(self, x: T, y: T, z: T) -> None:
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self) -> str:
        return f'Point({self.x}, {self.y}, {self.z})'

    def __add__(self, other: Point) -> Point:
        return Point(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: Point) -> Point:
        return Point(self.x - other.x, self.y - other.y, self.z - other.z)

class IVec(Point):
    def __init__(self, x: int, y: int, z: int):
        super().__init__(x, y, z)

class RVec(Point):
#class RVec(Point[float]):
    def __init__(self, x: float, y: float, z: float):
        super().__init__(x, y, z)

class Prop(object):
    """
    A class representing a prognostic variable, such as kinetic Energy,
    total Enery and temperature, in an MD simulation.

    Attributes
    ----------
    val : float
        the property value
    sum : float
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
        return (self.sum, self.sum2)

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
    x  : int
        the property value
    y  : int
        the sum of the property
    """

    def __init__(self, x : int=0, y : int=0, z : int=0):
        self.x = int(x)
        self.y = int(y)
        self.z = int(z)

    def __add__(self, other):
        """Return component vector addition.
        """
        return VecI((self.x - other.x), (self.y - other.y), (self.z - other.z))

    def truediv(a, b: VecR) -> VecR:
        "Same as a / b."
        return VecR(a.x / b.x, a.y / b.y, a.z / b.z)

    __truediv__ = truediv

    def mul(a: object, b: object) -> 'VecI':
        return VecI(a.x * b.x, a.y * b.y, a.z * b.z)

    rmul = mul # commutative operation

    def __repr__(self):
        return 'VecI({self.x}, {self.y}, {self.z})'.format(self=self)

    def vol(self) -> int:
        return self.x * self.y * self.z

    def vc_to_list_index(self, cells: 'VecI') -> int:
        """Translate vector cell index to scalar index using column-major order.
        """
        return (self.z * cells.y + self.y) * cells.x + self.x

class VecR(object):
    """
    A class representing a point in a 2-Dimensional Cartessian coordinate
    system of reference during an MD simulation.

    Attributes
    ==========
    x  : float
        the x coordinate
    y  : float
        the y coordinate
    z  : float
        the z coordinate
    """
    def __init__(self, x : float=0., y : float=0., z : float=0.):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, rhs):
        """Return the vector sum."""
        return VecR(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)

    __radd__ = __add__

    def __iadd__(self, other) -> VecR:
        return VecR(self.x + other.x, self.y + other.y, self.z + other.z)

    def truediv(a, b) -> VecR:
        "Same as a / b."
        return VecR(a.x / b.x, a.y / b.y, a.z / b.z)

    __truediv__ = truediv

    def __sub__(self, other):
        """Return relative vector difference."""
        return VecR(self.x - other.x, self.y - other.y, self.z - other.z)

    def __imul__(self, rhs: 'VecR') -> 'VecR':
        return VecR(self.x * rhs.x, self.y * rhs.y, self.z * rhs.z)

    def mul(a: object, b: object) -> 'VecR':
        return VecR(a.x * b.x, a.y * b.y, a.z * b.z)

    rmul = mul # commutative operation

    def __mul__(self, rhs: float) -> 'VecR':
        return VecR(self.x * rhs, self.y * rhs, self.z * rhs)

    __rmul__ = __mul__ # commutative operation

    def vcsum(self):
        return self.x + self.y + self.z

    def wrap(self, region):
        """
        Parameters
        ==========
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

        # Wrap the z-coordinate
        if self.z >= 0.5 * region.z:
            self.z -= region.z
        elif self.z < -0.5 * region.z:
            self.z += region.z

        return self

    def zero(self):
        """Zero the vector components.
        """
        self.x = 0.
        self.y = 0.
        self.z = 0.

        return self

    def __repr__(self):
        return 'VecR({self.x}, {self.y}, {self.z})'.format(self=self)

class Mol(object):
    """
    A class representing a molecule, i.e., position: r=(x,y), velocity: rv=(vx,vy),
    and acceleration: ra=(ax,ay), in an MD simulation.

    Attributes
    ----------
    r  : float
        the position vector
    rv  : float
        the velocity vector
    ra : float
        the acceleration vector

    Methods
    -------
    r_diff(other)
        Computes vector difference between self and other.
    r_wrap(region)
        Wraps the position vector to within boundary conditions.
    ra_zero()
        Zeros the acceleration vector.
    update_coordinates(integration_scheme, *args)
        Integrate the coordinates using integration scheme function.
    update_velocities(integration_scheme, *args)
        Integrate the velocities using integration scheme function.
    """
    def __init__(self, r: VecR=None, rv: VecR=None, ra: VecR=None):
        self.r  = r
        self.rv = rv
        self.ra = ra

        # Predictor-Corrector support
        self.ro  = None
        self.rvo = None
        self.ra1 = None
        self.ra2 = None

    def __repr__(self):
        return 'Mol({self.r}, {self.rv}, {self.ra})'.format(self=self)

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
        self.r.wrap(region)

        return self

    def ra_zero(self):
        self.ra.zero()

        return self

    def ra_sadd(self, s, v):
        """Scale molecular acceleration components.
        """
        self.ra.x += (s * v.x)
        self.ra.y += (s * v.y)
        self.ra.z += (s * v.z)

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
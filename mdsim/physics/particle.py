"""A Cartesian-coordinate point class.
"""
from __future__ import annotations
from typing import Generic, TypeAlias, TypeVar

T = TypeVar('T', int, float)
class Particle(Generic[T]):
    """Class representing a point in an inertial Cartesian coordinate system in
    which the components of ^r are (x, y, z)."""
    def __init__(self, x: T, y: T, z: T) -> None:
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self) -> str:
        return f'Particle({self.x}, {self.y}, {self.z})'

    def __add__(self, other: Particle) -> Particle:
        return Particle(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: Particle) -> Particle:
        return Particle(self.x - other.x, self.y - other.y, self.z - other.z)

def vec_sum[T: Particle](v1: T, v2: T) -> T:
    return v1 + v2
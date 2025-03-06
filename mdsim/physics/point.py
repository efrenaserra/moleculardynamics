"""A Cartesian-Coordinate point class.
"""
from __future__ import annotations
from typing import TypeVar

from ..core import Point

T = TypeVar('Point')

def sum[T](v1: T, v2: T) -> T:
    return Point(v1.x + v2.x, v1.x + v2.x, v1.z + v2.z)

def vec_sum[T](v1: T, v2: T) -> T:
    return v1 + v2
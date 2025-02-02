"""A Cartesian-Coordinate point class.
"""
class Point:
    def __init__(self, dimensions: list[int]):
        self.x =
        self.y =
        self.z =

    def coords(self):
        return tuple(self.x, self.y, self.z)
    
    def __sub__(self, other):
        return calc.distance(self, other)
    
    def __rsub__(self, other):
        self.__sub__(other)
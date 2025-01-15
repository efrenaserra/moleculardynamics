# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 21:32:12 2019

@author: Efren Antonio Serra (Read)
"""
class _Cells(object):
    def __init__(self, x: int=0, y: int=0):
        self.x = x
        self.y = y

class _Region(object):
    def __init__(self, x: float=0., y: float=0.):
        self.x = x
        self.y = y

class _VecI(object):
    def __init__(self, x: int=0, y: int=0):
        self.x = x
        self.y = y
    def __add__(self, rhs) -> '_VecI':
        return _VecI(self.x+rhs.x, self.y+rhs.y)

    def __repr__(self) -> str:
        return '_VecI({self.x}, {self.y})'.format(self=self)

    def vc_to_list_index(self, cells: _Cells) -> int:
        return (self.y * cells.x) + self.x

class _VecR(object):
    def __init__(self, x: float=0, y: float=0):
        self.x = x
        self.y = y
    def __add__(self, rhs) -> '_VecR':
        return _VecI(self.x+rhs.x, self.y+rhs.y)

    def __repr__(self) -> str:
        return '_VecR({self.x}, {self.y})'.format(self=self)

    def zero(self):
        self.x = self.y = 0.

vcOff = [_VecI(0,0), _VecI(1,0), _VecI(1,1), _VecI(0,1), _VecI(-1,1)]

def test_vc_offset():
    cells = _Cells(2,2)
    region = _Region(3.118,3.118)
    rshift = _VecR()

    for m1y in range(cells.y):
        for m1x in range(cells.x):
            m1v=_VecI(m1x,m1y)
            print("Cell vector",m1v," scalar cell index ",m1v.vc_to_list_index(cells))
            for offset in range(5):
                rshift.zero()
                m2v=m1v+vcOff[offset]
                print("Cell vector before BC",m2v," neighbor scalar cell index ",m2v.vc_to_list_index(cells))
                """Boundary conditions"""
                if m2v.x >= cells.x:
                    m2v.x = 0
                    rshift.x = region.x
                elif m2v.x < 0:
                    m2v.x = cells.x - 1
                    rshift.x = -region.x

                if m2v.y >= cells.y:
                    m2v.y = 0
                    rshift.y = region.y
                elif m2v.y < 0:
                    m2v.y = cells.y - 1
                    rshift.y = -region.y

                print("Cell vector after BC",m2v," neighbor scalar cell index ",m2v.vc_to_list_index(cells),rshift)

if __name__ == "__main__":
    test_vc_offset()
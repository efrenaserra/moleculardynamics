# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 13:53:43 2019

@author: Efren A. Serra
"""

import math

class Prop(object):
    def __init__(self):
        self.val : float = 0.
        self.sum : float = 0.
        self.sum2: float = 0.
    
    def accum(self):
        self.sum += self.val
        self.sum2 += math.sqrt(self.val)
    
    def zero(self):
        self.sum = self.sum2 = 0.

    def avg(self, n: int):
        self.sum /= n
        self.sum2 = math.sqrt(max([self.sum2 / n - math.sqrt(self.sum), 0.]))
    
    def est(self):
        return (self.sum, self.sum2)

    def __repr__(self):
        return "<val: %f; sum: %f; sum2: %f>"%(self.val,self.sum,self.sum2)

def test_prop():
    totEnergy = Prop()
    print(totEnergy)

if __name__ == "__main__":
    test_prop()
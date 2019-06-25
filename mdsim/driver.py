# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 14:12:58 2019

@author: Efren A. Serra
"""

def RunMDSim():
    GetNameList()
    PrintNameList()
    SetParams()
    SetupJob()
    moreCycles = True
    while moreCycles:
        SingleStep()
        if stepCount >= stepLimit:
            moreCycles = False

if __name__ == "__main__":
    RunMDSim()
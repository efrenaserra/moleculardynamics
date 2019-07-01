# -*- coding: utf-8 -*-
"""Created on Tue Jun 25 14:12:58 2019

@author: Efren A. Serra
"""

import math, numpy as np, re, sys

from _globals import NDIM, _mdsim_globals, _namelist_converter
from _types   import (
        Mol, Prop, VecR, r_diff_vectorized, rv_diff_vectorized,\
        ra_diff_vectorized
        )

def AccumProps(icode: int):
    """Accumulate thermodynamics properties.
    Parameters
    ----------
    icount : int,
    """
    if icode == 0:
        _mdsim_globals['totEnergy'].zero()
        _mdsim_globals['kinEnergy'].zero()
        _mdsim_globals['pressure'].zero()
    elif icode == 1:
        _mdsim_globals['totEnergy'].accum()
        _mdsim_globals['kinEnergy'].accum()
        _mdsim_globals['pressure'].accum()        
    elif icode == 2:
        stepAvg = _mdsim_globals['stepAvg']
        _mdsim_globals['totEnergy'].avg(stepAvg)
        _mdsim_globals['kinEnergy'].avg(stepAvg)
        _mdsim_globals['pressure'].avg(stepAvg)  

def AllocArrays():
    nMol = _mdsim_globals['nMol']
    # the molecules
    _mdsim_globals['mol'] = \
    np.array([Mol(VecR(), VecR(), VecR()) for i in range(nMol)], dtype=Mol)

def ApplyBoundaryCond():
    return

def ComputeForces():
    fcVal = rr = rrCut = rri = rri3 = 0.
    j1 = j2 = n = 0

    rrCut = math.sqrt(_mdsim_globals['rCut'])

    mol = _mdsim_globals['mol']
    nMol = _mdsim_globals['nMol']

    for n in range(nMol):
        mol[n].ra.x = 0.
        mol[n].ra.y = 0.

    _mdsim_globals['uSum'] = 0.
    _mdsim_globals['virSum'] = 0.

    dr = r_diff_vectorized(mol[1:], mol[:-1])
    drv = rv_diff_vectorized(mol[1:], mol[:-1])
    dra = ra_diff_vectorized(mol[1:], mol[:-1])

    dr.size, drv.size, dra.size

def EvalProps():
    return

def GetNameList(fd: str):
    with open(fd, 'r') as f:
        pattern = re.compile(r'initUcell')
        for line in f:
            m = pattern.match(line)
            if m:
                k = 'initUcell'
                line = line[len(k):]
                (nx,ny) = line.split()
                # matrix of molecular unit cells
                _mdsim_globals[k] = \
                _namelist_converter[k](nx,ny)
            else:
                (k,v) = line.split()
                _mdsim_globals[k] = \
                _namelist_converter[k](v)

def LeapfrogStep(count: int):
    return

def PrintNameList(fd: object):
    print(_mdsim_globals, file=fd)

def PrintSummary(fd: object):
    nMol = _mdsim_globals['nMol']
    print("%5d %8.4f %7.4f %s %s %s"%(\
          _mdsim_globals['stepCount'],\
          _mdsim_globals['timeNow'],  \
          _mdsim_globals['vSum'].vcsum() / nMol,     \
          _mdsim_globals['totEnergy'].est(),\
          _mdsim_globals['kinEnergy'].est(),\
          _mdsim_globals['pressure'].est()), \
          file=fd)

def InitCoords():
    return

def InitVels():
    _mdsim_globals['vSum'] = VecR(x=0.,y=0.)

def InitAccels():
    mol = _mdsim_globals['mol']
    nMol = _mdsim_globals['nMol']
    for n in range(nMol):
        mol[n].ra.x = 0.
        mol[n].ra.y = 0.

def SetupJob():
    """Setup global variables prior to simulation.
    """
    AllocArrays()
    _mdsim_globals['stepCount'] = 0
    InitCoords()
    InitVels()
    InitAccels()
    AccumProps(0)

def SetParams():
    density = _mdsim_globals['density']
    initUcell = _mdsim_globals['initUcell']

    _mdsim_globals['rCut'] = math.pow(2.,1./6.)
    _mdsim_globals['region'] = \
    VecR(x=1. / math.sqrt(density) * initUcell[0], \
         y=1. / math.sqrt(density) * initUcell[1])

    # the total number of molecules
    nMol = \
    _mdsim_globals['initUcell'][0] * _mdsim_globals['initUcell'][1]
    _mdsim_globals['nMol'] = nMol
    _mdsim_globals['velMag'] = \
    math.sqrt(NDIM * (1. - 1. / nMol) * _mdsim_globals['temperature'])

    # initialize totEnery and kinEnergy properties
    _mdsim_globals['kinEnergy'] = Prop()
    _mdsim_globals['pressure']  = Prop()
    _mdsim_globals['totEnergy'] = Prop()

def SingleStep():
    _mdsim_globals['stepCount'] += 1
    _mdsim_globals['timeNow'] = \
    _mdsim_globals['stepCount'] * _mdsim_globals['deltaT']
    LeapfrogStep(1)
    ApplyBoundaryCond()
    ComputeForces()
    LeapfrogStep(2)
    EvalProps()
    AccumProps(1)
    if _mdsim_globals['stepCount'] % _mdsim_globals['stepAvg'] == 0:
        AccumProps(2)
        PrintSummary(sys.stdout)
        AccumProps(0)

def RunMDSim(argv: list):
    GetNameList(argv[1])
    PrintNameList(sys.stdout)
    SetParams()
    SetupJob()
    moreCycles = True
    while moreCycles:
        SingleStep()
        if _mdsim_globals['stepCount'] >= _mdsim_globals['stepLimit']:
            moreCycles = False

if __name__ == "__main__":
    RunMDSim(['driver.py', 'pr_02_01.in'])
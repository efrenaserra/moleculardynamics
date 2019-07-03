# -*- coding: utf-8 -*-
"""Created on Tue Jun 25 14:12:58 2019

@author: Efren A. Serra
"""

import math, numpy as np, re, sys

from _globals import NDIM, _mdsim_globals, _namelist_converter
from _types   import (
        Mol, Prop, VecR
        )
from _vfunctions import (
        r_diff, \
        r_diff_vectorized, \
        r_wrap, \
        r_wrap_vectorized, \
        rv_add_vectorized, \
        rv_diff_vectorized, \
        rv_dot_vectorized, \
        rv_rand_vectorized, \
        rv_scale_vectorized, \
        rv_sadd_vectorized, \
        ra_diff_vectorized, \
        ra_zero_vectorized, \
        leapfrog_integrate_rv_vectorized, \
        leapfrog_integrate_r_vectorized,
        vecr_wrap, \
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
    """
    """
    nMol = _mdsim_globals['nMol']
    # the molecules
    mol = \
    np.array([Mol(VecR(), VecR(), VecR()) for i in range(nMol)], dtype=Mol)
    _mdsim_globals['mol'] = mol[np.newaxis, :]

def ApplyBoundaryCond():
    """
    """
    mol    = _mdsim_globals['mol']
    region = _mdsim_globals['region']

    r_wrap_vectorized(mol, region)

def ComputeForces():
    """Compute the MD forces by evaluating the LJ potential
    """
    fcVal : float = 0.
    rr : float    = 0.
    rrCut : float = 0.
    rri : float   = 0.
    rri3 : float  = 0.
    j1 : int = 0
    j2 : int = 0

    mol    = _mdsim_globals['mol']
    nMol   = _mdsim_globals['nMol']
    region = _mdsim_globals['region']
    rrCut  = math.sqrt(_mdsim_globals['rCut'])

    ra_zero_vectorized(mol)
    _mdsim_globals['uSum'] = 0.
    _mdsim_globals['virSum'] = 0.

    # compute cross differences
    for j1, a in enumerate(mol[0, :nMol-1]):
        for j2, b in enumerate(mol[0, j1+1:]):
            dr = r_diff(a, b)
            dr = vecr_wrap(dr, region)

def EvalProps():
    """Evaluate thermodynamic properties
    """
    density  = _mdsim_globals['density']
    mol  = _mdsim_globals['mol']
    nMol = _mdsim_globals['nMol']
    uSum = _mdsim_globals['uSum']
    vSum = _mdsim_globals['vSum']
    virSum = _mdsim_globals['virSum']

    vSum.zero()
    vvSum : float = 0.
    rv_add_vectorized(vSum, mol)
    vvSum = rv_dot_vectorized(mol, mol)
    
    _mdsim_globals['kinEnergy'].val = 0.5 * vvSum / nMol
    _mdsim_globals['totEnergy'].val = \
    _mdsim_globals['kinEnergy'].val + uSum / nMol
    _mdsim_globals['pressure'].val = density * (vvSum + virSum) / (nMol * NDIM)

def GetNameList(fd: str):
    """
    Parameters
    ----------
    fd : str, the filename
    """
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

def LeapfrogStep(part: int):
    """
    Parameters
    ----------
    part : int, 
    """
    deltaT : float = _mdsim_globals['deltaT']
    mol            = _mdsim_globals['mol']

    if part == 1:
        # integrate velocities
        leapfrog_integrate_rv_vectorized(mol, 0.5 * deltaT)
        # integrate coordinates
        leapfrog_integrate_r_vectorized(mol, deltaT)
    else:
        # integrate velocities
        leapfrog_integrate_rv_vectorized(mol, 0.5 * deltaT)

def PrintNameList(fd: object):
    """
    Parameters
    ----------
    fd : object, 
    """
    print(_mdsim_globals, file=fd)

def PrintSummary(fd: object):
    """
    Parameters
    ----------
    fd : object, 
    """
    nMol = _mdsim_globals['nMol']
    print("%5d %8.4f %7.4f %s %s %s"%(           \
          _mdsim_globals['stepCount'],           \
          _mdsim_globals['timeNow'],             \
          _mdsim_globals['vSum'].vcsum() / nMol, \
          _mdsim_globals['totEnergy'].est(),     \
          _mdsim_globals['kinEnergy'].est(),     \
          _mdsim_globals['pressure'].est()),     \
          file=fd)

def InitCoords():
    return

def InitVels():
    mol    = _mdsim_globals['mol']
    nMol   = _mdsim_globals['nMol']
    velMag = _mdsim_globals['velMag']

    vSum = VecR(x=0., y=0.)
    rv_rand_vectorized(mol)
    rv_scale_vectorized(mol, velMag)
    rv_add_vectorized(vSum, mol)
    _mdsim_globals['vSum'] = vSum
    # scale molecular velocities
    rv_sadd_vectorized(mol, - 1. / nMol, vSum)

def InitAccels():
    mol = _mdsim_globals['mol']
    ra_zero_vectorized(mol)

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
    nMol = initUcell[0] * initUcell[1]
    _mdsim_globals['nMol'] = nMol
    _mdsim_globals['velMag'] = \
    math.sqrt(NDIM * (1. - 1. / nMol) * _mdsim_globals['temperature'])

    # initialize kinEnery, pressure and totEnergy properties
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
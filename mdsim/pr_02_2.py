# -*- coding: utf-8 -*-
"""Created on Tue Jun 25 14:12:58 2019

/* [[pr_02_2 - velocity distribution]] */

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

import math, numpy as np, re, sys
import matplotlib.pyplot as plt

from _globals import NDIM, _mdsim_globals, _namelist_converter
from _types   import (
        Mol, Prop, VecR
        )

from _leapfrog_functions import (
        leapfrog_update_coordinates,\
        leapfrog_update_velocities, \
        )

from _vec_functions import (
        rv_add,   \
        rv_dot,   \
        rv_rand,  \
        rv_sadd,  \
        rv_scale, \
        vecr_div, \
        vecr_dot, \
        vecr_mul, \
        vecr_sadd,\
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
    """Allocate array of molecules.
    """
    nMol = _mdsim_globals['nMol']
    sizeHistVel = _mdsim_globals['sizeHistVel']

    # The molecules
    mol = \
    np.array([Mol() for i in range(nMol)], dtype=Mol)
    _mdsim_globals['mol'] = mol

    # The velocity histogram
    histVel = \
    np.array([0.0 for i in range(sizeHistVel)], dtype=float)
    _mdsim_globals['histVel'] = histVel

def ApplyBoundaryCond():
    """Apply periodic boundary conditions.
    """
    mol    = _mdsim_globals['mol']
    region = _mdsim_globals['region']

    for i, m in enumerate(mol):
        m.r_wrap(region)

def ComputeForces():
    """Compute the MD forces by evaluating the LJ potential
       using all-pairs interactions.
    """
    j1 : int = 0
    j2 : int = 0
    fcVal : float = 0.
    rr : float = 0.
    rrCut : float = 0.
    rri : float   = 0.
    rri3 : float  = 0.

    mol    = _mdsim_globals['mol']
    nMol   = _mdsim_globals['nMol']
    region = _mdsim_globals['region']
    rrCut  = _mdsim_globals['rCut'] * _mdsim_globals['rCut']

    for m in mol:
        m.ra_zero()

    _mdsim_globals['uSum'] = 0.
    _mdsim_globals['virSum'] = 0.

    # Compute all-pairs interactions
    _mdsim_globals['uSum'] = 0.
    for j1 in range(nMol-1):
        a = mol[j1]
        for j2 in range(j1+1, nMol):
            b = mol[j2]
            dr = a.r_diff(b)
            dr.wrap(region)
            rr = vecr_dot(dr, dr)
            if rr < rrCut:
                rri = 1. / rr
                rri3 = rri ** 3
                fcVal = 48.0 * rri3 * (rri3 - 0.5) * rri

                # Molecule at: j1
                a.ra_sadd(fcVal, dr)
                # Molecule at: j2
                b.ra_sadd(-fcVal, dr)

                _mdsim_globals['uSum'] += 4. * rri3 * (rri3 - 1.) + 1.
                _mdsim_globals['virSum'] += fcVal * rr

def DisplayTrajectories():
    mol   = _mdsim_globals['mol']
    for n, m in enumerate(mol):
        if n % 25 == 0:
            plt.scatter(m.r.x, m.r.y, marker='o', s=5, c='red')
        else:
            plt.scatter(m.r.x, m.r.y, marker='o', s=5, c='black')
    plt.show()

def EvalProps():
    """Evaluate thermodynamic properties
    """
    density= _mdsim_globals['density']
    mol    = _mdsim_globals['mol']
    nMol   = _mdsim_globals['nMol']
    uSum   = _mdsim_globals['uSum']
    virSum = _mdsim_globals['virSum']

    vSum = VecR()
    vvSum : float = 0.
    for m in mol:
        rv_add(vSum, m)
        vvSum += rv_dot(m, m)

    _mdsim_globals['vSum'] = vSum
    _mdsim_globals['kinEnergy'].val = 0.5 * vvSum / nMol
    _mdsim_globals['totEnergy'].val = \
    _mdsim_globals['kinEnergy'].val + uSum / nMol
    _mdsim_globals['pressure'].val = density * (vvSum + virSum) / (nMol * NDIM)

def EvalVelDist():
    j : int = 0
    deltaV : float = 0.0
    histSum : float = 0.0

    histVel     = _mdsim_globals['histVel']
    mol         = _mdsim_globals['mol']
    rangeVel    = _mdsim_globals['rangeVel']
    sizeHistVel = _mdsim_globals['sizeHistVel']

    deltaV = rangeVel / sizeHistVel
    for m in mol:
        j =  int(math.sqrt(rv_dot(m,m)) / deltaV)
        histVel[min([j, sizeHistVel - 1])] += 1

    _mdsim_globals['countVel'] += 1
    if _mdsim_globals['countVel'] == _mdsim_globals['limitVel']:
        histSum = np.sum(histVel)
        histVel /= histSum
        _mdsim_globals['hFucnction'] = 0.0
        for j in range(sizeHistVel):
            if histVel[j] > 0.0:
                _mdsim_globals['hFunction'] += histVel[j] * math.log(histVel[j] / ((j+0.5) * deltaV))
        PrintVelDist(sys.stdout)
        _mdsim_globals['countVel'] = 0

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
                # Matrix of molecular unit cells
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
        for i, m in enumerate(mol):
            # Integrate velocities
            m.update_velocities(leapfrog_update_velocities, 0.5 * deltaT)

            # Integrate coordinates
            m.update_coordinates(leapfrog_update_coordinates, deltaT)
    else:
        for i, m in enumerate(mol):
            # Integrate velocities
            m.update_velocities(leapfrog_update_velocities, 0.5 * deltaT)

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
    totEnergy='%7.4f %7.4f'%_mdsim_globals['totEnergy'].est()
    kinEnergy='%7.4f %7.4f'%_mdsim_globals['kinEnergy'].est()
    pressure='%7.4f %7.4f'%_mdsim_globals['pressure'].est()
    print("%5d %8.4f %7.4f %s %s %s"%(\
          _mdsim_globals['stepCount'],\
          _mdsim_globals['timeNow'],  \
          _mdsim_globals['vSum'].vcsum() / nMol,\
          totEnergy, \
          kinEnergy, \
          pressure),  \
          file=fd)

def PrintVelDist(fd: object):
    """
    Parameters
    ----------
    fd : object, 
    """
    histVel     = _mdsim_globals['histVel']
    rangeVel    = _mdsim_globals['rangeVel']
    sizeHistVel = _mdsim_globals['sizeHistVel']

    print("vdist (%.3f)"%(_mdsim_globals['timeNow']), file=fd)
    bins = []
    for n in range(sizeHistVel):
        vBin = (n + 0.5) * rangeVel / sizeHistVel
        bins.append(vBin)
        print("%8.3f %8.3f"%(vBin, histVel[n]), file=fd)
    print("hfun: (%8.3f %8.3f)"%(_mdsim_globals['timeNow'],_mdsim_globals['hFunction']), file=fd)
    plt.scatter(bins, histVel, marker='+')
    plt.ylim(0.0)
    plt.ylabel('f(|v|)')
    plt.xlabel('|v|')
    plt.show()

def InitCoords():
    """Initialize the molecular coordinates
    """
    mol       = _mdsim_globals['mol']
    region    = _mdsim_globals['region']
    initUcell = _mdsim_globals['initUcell']

    gap = vecr_div(region, initUcell)
    for nx in range(initUcell.x):
        for ny in range(initUcell.y):
            c = VecR(x=nx + 0.5, y=ny + 0.5)
            c = vecr_mul(c, gap)
            c = vecr_sadd(c, -0.5, region)
            mol[nx * initUcell.x + ny].r = c

def InitVels():
    """Initialize the molecular velocities
    """
    mol    = _mdsim_globals['mol']
    nMol   = _mdsim_globals['nMol']
    velMag = _mdsim_globals['velMag']

    vSum = VecR(x=0., y=0.)
    for m in mol:
        m.rv = VecR()
        rv_rand(m)
        rv_scale(m, velMag)
        rv_add(vSum, m)

    _mdsim_globals['vSum'] = vSum
    # Scale molecular velocities
    for m in mol:
        rv_sadd(m, -1. / nMol, vSum)

def InitAccels():
    """Initialize the molecular accelerations
    """
    mol = _mdsim_globals['mol']
    for m in mol:
        m.ra = VecR()
        m.ra_zero()

def SetupJob():
    """Setup global variables prior to simulation.
    """
    AllocArrays()
    _mdsim_globals['stepCount'] = 0
    InitCoords()
    InitVels()
    InitAccels()
    AccumProps(0)
    _mdsim_globals['countVel'] = 0
    _mdsim_globals['hFunction'] = 0.0

def SetParams():
    density = _mdsim_globals['density']
    initUcell = _mdsim_globals['initUcell']

    _mdsim_globals['rCut'] = math.pow(2.,1./6.)
    _mdsim_globals['region'] = \
    VecR(x=1. / math.sqrt(density) * initUcell.x, \
         y=1. / math.sqrt(density) * initUcell.y)

    # The total number of molecules
    nMol = initUcell.x * initUcell.y
    _mdsim_globals['nMol'] = nMol
    _mdsim_globals['velMag'] = \
    math.sqrt(NDIM * (1. - 1. / nMol) * _mdsim_globals['temperature'])

    # Initialize kinEnery, pressure and totEnergy properties
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
    if _mdsim_globals['stepCount'] >= _mdsim_globals['stepEquil'] and \
    (_mdsim_globals['stepCount'] - _mdsim_globals['stepEquil']) % _mdsim_globals['stepVel'] == 0:
        EvalVelDist()
    if _mdsim_globals['stepCount'] % _mdsim_globals['stepAvg'] == 0:
        AccumProps(2)
        PrintSummary(sys.stdout)
        DisplayTrajectories()
        AccumProps(0)

def RunMDSim(argv: list):
    GetNameList(argv[0])
    PrintNameList(sys.stdout)
    SetParams()
    SetupJob()
    moreCycles = True
    while moreCycles:
        SingleStep()
        if _mdsim_globals['stepCount'] >= _mdsim_globals['stepLimit']:
            moreCycles = False

if __name__ == "__main__":
    RunMDSim(['pr_02_2.in'])
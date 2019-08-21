# -*- coding: utf-8 -*-
"""Created on Tue Jun 25 14:12:58 2019

/* [[pr_03_2 - neighbor list and leapfrog]] */

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

NDIM : int = 3
N_OFFSET : int = 14
import math, numpy as np, re, sys

from _globals import _mdsim_globals, _namelist_converter
from _types   import (
        Mol, Prop, VecI, VecR
        )

from _leapfrog_functions import (
        leapfrog_update_coordinates,
        leapfrog_update_velocities,
        )

from _3d_rand import (
        rv_rand
        )

from _vec_functions import (
        rv_add,
        rv_dot,
        rv_sadd,
        rv_scale,
        vecr_dot,
        _VCell_wrap_all,
        )

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

def SingleStep():
    _mdsim_globals['stepCount'] += 1
    _mdsim_globals['timeNow'] = \
    _mdsim_globals['stepCount'] * _mdsim_globals['deltaT']

    LeapfrogStep(1)
    ApplyBoundaryCond()
    if _mdsim_globals['nebrNow']:
        _mdsim_globals['nebrNow'] = False
        _mdsim_globals['dispHi'] = 0.
        BuildNebrList()
    ComputeForces()
    LeapfrogStep(2)
    EvalProps()
    if _mdsim_globals['stepCount'] < _mdsim_globals['stepEquil']:
        AdjustInitTemp()
    AccumProps(1)
    if _mdsim_globals['stepCount'] % _mdsim_globals['stepAvg'] == 0:
        AccumProps(2)
        PrintSummary(sys.stdout)
        AccumProps(0)

def SetupJob():
    """Setup global variables prior to simulation.
    """
    AllocArrays()
    _mdsim_globals['stepCount'] = 0
    InitCoords()
    InitVels()
    InitAccels()
    AccumProps(0)
    _mdsim_globals['kinEnInitSum'] = 0.0
    _mdsim_globals['nebrNow'] = True

def SetParams():
    density    = _mdsim_globals['density']
    initUcell  = _mdsim_globals['initUcell']
    nebrTabFac = _mdsim_globals['nebrTabFac']
    rNebrShell = _mdsim_globals['rNebrShell']

    rCut = math.pow(2., 1./6.)
    _mdsim_globals['rCut'] = rCut
    region =\
    VecR(x=1. / math.pow(density, 1./3.) * initUcell.x, \
         y=1. / math.pow(density, 1./3.) * initUcell.y, \
         z=1. / math.pow(density, 1./3.) * initUcell.z)
    _mdsim_globals['region'] = region

    # The total number of molecules
    nMol = initUcell.x * initUcell.y * initUcell.z
    _mdsim_globals['nMol'] = nMol

    _mdsim_globals['velMag'] = \
    math.sqrt(NDIM * (1. - 1. / nMol) * _mdsim_globals['temperature'])

    _mdsim_globals['cells'] = \
    VecI(x=(1. / (rCut + rNebrShell)) * region.x, \
         y=(1. / (rCut + rNebrShell)) * region.y, \
         z=(1. / (rCut + rNebrShell)) * region.z)

    # The neighbor list size
    _mdsim_globals['nebrTabMax'] = nebrTabFac * nMol 
    
    # Initialize kinEnery, pressure and totEnergy properties
    _mdsim_globals['kinEnergy'] = Prop()
    _mdsim_globals['pressure']  = Prop()
    _mdsim_globals['totEnergy'] = Prop()

def AllocArrays():
    """Allocate array of molecules.
    """
    nMol       = _mdsim_globals['nMol']
    nebrTabMax = _mdsim_globals['nebrTabMax']

    # The molecules
    _mdsim_globals['mol'] = \
    np.array([Mol() for i in range(nMol)], dtype=Mol)

    # The cells
    cells : VecI = _mdsim_globals['cells']
    nCells = cells.vol()
    _mdsim_globals['cellList'] = \
    np.array([-1 for i in range(nCells + nMol)], dtype=int)

    # The neighbor table
    _mdsim_globals['nebrTab'] = \
    np.array([-1 for i in range(2 * nebrTabMax)], dtype=int)

def BuildNebrList():
    cc : VecI =  None
    m1v: VecI =  None
    m2v: VecI =  None
    # Mapping from cell vector to list index:
    #     (c.z * cells.y + c.y) * cells.x + c.x -> 0
    #
    # With cells = VecI(3, 3, 3), list indices are:
    vOff = [
            VecI(0,0,0),   # (0 * 3 + 0) * 3 + 0 -> 0
            VecI(1,0,0),   # (0 * 3 + 0) * 3 + 1 -> 1
            VecI(1,1,0),   # (0 * 3 + 1) * 3 + 1 -> 4
            VecI(0,1,0),   # (0 * 3 + 1) * 3 + 0 -> 3
            VecI(-1,1,0),  # (0 * 3 + 1) * 3 - 1 -> 2
            VecI(0,0,1),   # (1 * 3 + 0) * 3 + 0 -> 9
            VecI(1,0,1),   # (1 * 3 + 0) * 3 + 1 -> 10
            VecI(1,1,1),   # (1 * 3 + 1) * 3 + 1 -> 13
            VecI(0,1,1),   # (1 * 3 + 1) * 3 + 0 -> 12
            VecI(-1,1,1),  # (1 * 3 + 1) * 3 - 1 -> 11
            VecI(-1,0,1),  # (1 * 3 + 0) * 3 - 1 -> 8
            VecI(-1,-1,1), # (1 * 3 - 1) * 3 - 1 -> 5
            VecI(0,-1,1),  # (1 * 3 - 1) * 3 + 0 -> 6
            VecI(1,-1,1),  # (1 * 3 - 1) * 3 + 1 -> 7
            ]
    c  : int = 0
    j1 : int = 0
    j2 : int = 0

    cells    = _mdsim_globals['cells']
    cellList = _mdsim_globals['cellList']
    nCells   = cells.vol()

    mol      = _mdsim_globals['mol']
    nMol     = _mdsim_globals['nMol']
    rCut     = _mdsim_globals['rCut']
    region   = _mdsim_globals['region']
    rNebrShell = _mdsim_globals['rNebrShell']
    nebrTab  = _mdsim_globals['nebrTab']
    nebrTabMax = _mdsim_globals['nebrTabMax']
    rrNebr   =  (rCut + rNebrShell) * (rCut + rNebrShell)

    invWid: VecR = cells / region
    for n in range(nMol, nMol + nCells):
        cellList[n] = -1

    for n in range(nMol):
        rs = mol[n].r + (0.5 * region)
        cc = VecI(rs.x * invWid.x, rs.y * invWid.y, rs.z * invWid.z)
        c  = cc.vc_to_list_index(cells) + nMol
        cellList[n] = cellList[c]
        cellList[c] = n

    _mdsim_globals['nebrTabLen'] = 0
    rshift = VecR()
    for m1z in range(cells.z):
        for m1y in range(cells.y):
            for m1x in range(cells.x):
                """Vector cell index to which this molecule belongs."""
                m1v = VecI(m1x,m1y,m1z)
                """Translate vector cell index, m1v, to cell list index."""
                m1 = m1v.vc_to_list_index(cells) + nMol
                for offset in range(N_OFFSET):
                    """Vector cell relative to m1v."""
                    m2v = m1v + vOff[offset]
                    """ Periodic boundary-conditions/shift coordinates."""
                    rshift.zero()
                    _VCell_wrap_all(m2v, cells, rshift, region)
                    """Translate vector cell index, m2v, to cell list index."""
                    m2 = m2v.vc_to_list_index(cells) + nMol
                    j1 = cellList[m1]
                    while j1 >= 0:
                        j2 = cellList[m2]
                        while j2 >= 0:
                            if m1 != m2 or j2 < j1:
                                a = mol[j1]
                                b = mol[j2]
                                dr = a.r_diff(b)
                                dr -= rshift
                                rr = vecr_dot(dr, dr)
                                if rr < rrNebr:
                                    if _mdsim_globals['nebrTabLen'] >= nebrTabMax:
                                        sys.exit(1)
                                    nebrTabLen = _mdsim_globals['nebrTabLen']
                                    nebrTab[2 * nebrTabLen] = j1
                                    nebrTab[2 * nebrTabLen + 1] = j2
                                    _mdsim_globals['nebrTabLen'] += 1
                            j2 = cellList[j2]
                        j1 = cellList[j1]

def ComputeForces():
    """Compute the MD forces by evaluating the LJ potential
       using cells subdivisions. The simulation area, i.e., region,
       is divided into a lattice of small cells and the cell edges all exceed
       rcut. For this implementation:
           rca = La/Lca, where Lca = FLOOR(La/rcut) (a = x, y, z), and La is
       simulation box length in the a direction (region).

           +---------------------------------------------+
           |          |                       |          |
           |       -1 | 0                Lcx-1|Lcx       |
           |      +---+---+---+---+---+---+---+---+      |
           |      |   |                       |   |      |
           +------+---+---+---+---+---+---+---+---+------+

    """
    j1    : int = 0
    j2    : int = 0

    fcVal : float = 0.
    uVal  : float = 0.
    rr    : float = 0.
    rrCut : float = 0.
    rri   : float = 0.
    rri3  : float = 0.

    mol        = _mdsim_globals['mol']
    nebrTab    = _mdsim_globals['nebrTab']
    nebrTabLen = _mdsim_globals['nebrTabLen']
    region     = _mdsim_globals['region']
    rrCut      = _mdsim_globals['rCut'] * _mdsim_globals['rCut']

    for m in mol:
        m.ra_zero()

    _mdsim_globals['uSum']   = 0.
    _mdsim_globals['virSum'] = 0.

    # Compute cells interactions
    for n in range(nebrTabLen):
        j1 = nebrTab[2 * n]
        j2 = nebrTab[2 * n + 1]
        a = mol[j1]
        b = mol[j2]
        dr = a.r_diff(b)
        dr.wrap(region)
        rr = vecr_dot(dr, dr)
        if rr < rrCut:
            rri = 1. / rr
            rri3 = rri ** 3
            fcVal = 48.0 * rri3 * (rri3 - 0.5) * rri
            uVal = 4. * rri3 * (rri3 - 1.) + 1.
            # Molecule at: j1
            a.ra_sadd(fcVal, dr)
            # Molecule at: j2
            b.ra_sadd(-fcVal, dr)
            _mdsim_globals['uSum'] += uVal
            _mdsim_globals['virSum'] += fcVal * rr

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

def ApplyBoundaryCond():
    """Apply periodic boundary conditions.
    """
    mol    = _mdsim_globals['mol']
    region = _mdsim_globals['region']

    for m in mol:
        m.r_wrap(region)

def AdjustInitTemp():
    mol = _mdsim_globals['mol']
    _mdsim_globals['kinEnInitSum'] += _mdsim_globals['kinEnergy'].val

    if _mdsim_globals['stepCount'] % _mdsim_globals['stepInitlzTemp'] == 0:
        _mdsim_globals['kinEnInitSum'] /= _mdsim_globals['stepInitlzTemp']
        vFac = \
        _mdsim_globals['velMag'] / math.sqrt(2. * _mdsim_globals['kinEnInitSum'])
        for m in mol:
            rv_scale(m, vFac)
        _mdsim_globals['kinEnInitSum'] = 0.0

def InitCoords():
    """Initialize the molecular coordinates
    """
    mol       = _mdsim_globals['mol']
    region    = _mdsim_globals['region']
    initUcell = _mdsim_globals['initUcell']

#    gap = vecr_div(region, initUcell)
    gap = region / initUcell
    n = 0
    for nz in range(initUcell.z):
        for ny in range(initUcell.y):
            for nx in range(initUcell.x):
                c = VecR(x=nx + 0.5, y=ny + 0.5, z=nz + 0.5)
                c *= gap
                c += (-0.5 * region)
                mol[n].r = c
                n += 1

def InitVels():
    """Initialize the molecular velocities
    """
    mol    = _mdsim_globals['mol']
    nMol   = _mdsim_globals['nMol']
    velMag = _mdsim_globals['velMag']

    vSum = VecR(x=0., y=0., z=0.)
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

def EvalProps():
    """Evaluate thermodynamic properties
    """
    density = _mdsim_globals['density']
    deltaT  = _mdsim_globals['deltaT']
    mol     = _mdsim_globals['mol']
    nMol    = _mdsim_globals['nMol']
    uSum    = _mdsim_globals['uSum']
    virSum  = _mdsim_globals['virSum']

    vSum = VecR()
    vvMax : float = 0.
    vvSum : float = 0.
    for m in mol:
        rv_add(vSum, m)
        vv = rv_dot(m, m)
        vvSum += vv
        vvMax = max([vvMax, vv])

    _mdsim_globals['dispHi'] += math.sqrt(vvMax) * deltaT
    if _mdsim_globals['dispHi'] > 0.5 * _mdsim_globals['rNebrShell']:
        _mdsim_globals['nebrNow'] = True

    _mdsim_globals['vSum'] = vSum
    _mdsim_globals['kinEnergy'].val = 0.5 * vvSum / nMol
    _mdsim_globals['totEnergy'].val = \
    _mdsim_globals['kinEnergy'].val + uSum / nMol
    _mdsim_globals['pressure'].val = density * (vvSum + virSum) / (nMol * NDIM)

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
                (nx,ny,nz) = line.split()
                # Matrix of molecular unit cells
                _mdsim_globals[k] = \
                _namelist_converter[k](nx,ny,nz)
            else:
                (k,v) = line.split()
                _mdsim_globals[k] = \
                _namelist_converter[k](v)

def PrintNameList(fd: object):
    """
    Parameters
    ----------
    fd : object, the output file descriptor
    """
    print(_mdsim_globals, file=fd)

if __name__ == "__main__":
    RunMDSim(['pr_03_2.in'])
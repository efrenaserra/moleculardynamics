# -*- coding: utf-8 -*-
"""
Module defining global variables.

Created on Tue Jun 25 14:31:17 2019

@author: Efren A. Serra
"""
# Disallow reloading this module so as to prevent re-initializing these global
# variables
if '_is_loaded' in globals():
    raise RuntimeError('Reloading mdsim._globals is not allowed')
_is_loaded = True


"""Global variables
"""
NDIM : int = 2
_mdsim_globals = {}
_namelist_converter = {
    'deltaT'      : lambda x: float(x),
    'density'     : lambda x: float(x),
    'initUcell'   : lambda x,y: [int(x),int(y)],
    'limitVel'    : lambda x: int(x),
    'nebrTabFac'  : lambda x: int(x),
    'randSeed'    : lambda x: int(x),
    'rangeVel'    : lambda x: float(x),
    'rNebrShell'  : lambda x: float(x),
    'sizeHistVel' : lambda x: int(x),
    'stepAvg'     : lambda x: int(x),
    'stepEquil'   : lambda x: int(x),
    'stepLimit'   : lambda x: int(x),
    'stepVel'     : lambda x: int(x),
    'temperature' : lambda x: float(x),
    }
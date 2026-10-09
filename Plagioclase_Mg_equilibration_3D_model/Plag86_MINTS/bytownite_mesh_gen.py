#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:48:41 2019

@author: ejfm2
"""

# Script to generate plagioclase and anorthite mesh 

import pygmsh
import numpy as np
import meshio
import pandas as pd
import random
import sys

import xtal_geometries2 as xtal
import xtal_geometries_feldspars as fld

z = float(sys.argv[1]) # half length in microns
f_mesh = sys.argv[2]
f_res = float(sys.argv[3])

# Generate mesh resolution based on z axis

res = z/f_res


geom = pygmsh.built_in.Geometry()

#fsp = xtal.anorthite_743(geom, z, res)
fsp = fld.bytownite_dhz_fig210e(geom, z, res)

fsp_pv = fsp.phys_volume()

code = geom.get_code()

with open(f_mesh + '.geo','w') as out:
    out.write(code)

out.close()



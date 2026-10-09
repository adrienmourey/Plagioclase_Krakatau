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

z = float(sys.argv[1]) # half length in microns
f_mesh1 = sys.argv[2]
f_mesh2 = sys.argv[3]
f_mesh3 = sys.argv[4]
f_res = sys.argv[5]

#geom1 = pygmsh.built_in.Geometry()
geom2 = pygmsh.built_in.Geometry()
#geom3 = pygmsh.built_in.Geometry()

# Need to generate coordinate reference frame to position melt inclusion

# define plag size - scale along z axis (positive direction) - in microns 



#fsp1 = xtal.feldspar_656(geom1, z, float(f_res))
#fsp_pv = fsp1.phys_volume()

#code1 = geom1.get_code()

#with open(f_mesh1 + '.geo','w') as out1:
    #out1.write(code1)

#out1.close()


an1 = xtal.anorthite_743(geom2, z, float(f_res))
an_pv = an1.phys_volume()

code2 = geom2.get_code()

with open(f_mesh2 + '.geo','w') as out2:
    out2.write(code2)

out2.close() 


#fsp2 = xtal.feldspar_657(geom3, z, float(f_res))
#fsp_pv2 = fsp2.phys_volume()

#code3 = geom3.get_code()

#with open(f_mesh3 + '.geo','w') as out3:
    #out3.write(code3)

#out3.close() 


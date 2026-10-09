#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 11:46:25 2018

@author: ejfm2
"""

# Test script creating 3D spherical mesh of designated radius
# Extract central coordinate value for subsequent processing. 

from dolfin import *
import mshr
import numpy as np
import pandas as pd
#import KC_fO2 as kc
import math as m
import sys

set_log_active(False)

class Outer(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

def D(T_i, Xan, lnD0, cXan, cT_i):
    
    
    lnD = lnD0 + cXan*Xan + cT_i*T_i
    
    D = exp(lnD)*Constant(1e12)
    
    return D
    
    
def D1(T_i, Xan, lnD0, cXan, cT_i):
    
    
    lnD = lnD0 + cXan*Xan + cT_i*T_i
    
    D = np.exp(lnD)*1e12
    
    return D    
    
     
def mod_diff(t, T, Xan, A, rim_comp, core_comp, element):
    
    t *= (86400.0*365.0) # convert time into seconds # may have to start from years rather than days
    dt = t/nts # only 2000 time steps
    dT = Constant(dt)
    #P *= 1.0e8
    T += 273.0
    T_i = 1/T
    
    # Determine diffusion coefficient values based on element
    
    if element == 'Mg':
        D_x = np.array([-5.45, -7.983, -3.54e4]) # Van Orman
    if element == 'Ba': 
        D_x = np.array([-12.32, -3.287, -3.995e4])
    if element == 'Sr':
        D_x = np.array([-12.81, -5.712, -3.239e4])

    R = Constant(0.008314)
    
    Q = FunctionSpace(mesh, "CG", 1)
    
    
    # Define anorthite content as being constant
    An = Function(Q)
    
        
    # Construct weak form
    C0 = Function(Q)       # Composition at current time step
    C1 = TrialFunction(Q)  # Composition at next time step
    S = TestFunction(Q)
        
    theta = Constant(0.5)
    C_mid = theta*C1 + Constant((1.0-theta))*C0

    T = Constant(T)
    T_i = Constant(T_i)
    Xan = Constant(Xan)
    #aSiO2 = Constant(aSiO2)
    A = Constant(A)
    
    An.interpolate(Xan)
    
    # Plagioclase weak form
    
    D_pl = D(T_i, An, Constant(D_x[0]), Constant(D_x[1]), Constant(D_x[2]))
    
    F1 = (D_pl)/(R*T)*(A)*grad(An)
    F = S*(C1-C0)*dx + dT*(inner(grad(S), D_pl*grad(C_mid) - F1*(C_mid)))*dx
    
    #F = S*(C1-C0)*dx + dT*(inner(grad(S), D_pl*grad(C_mid)))*dx
            
    a = lhs(F)
    L = rhs(F)

    u0 = Constant(rim_comp)
    OB = Outer()

    Cbcs = DirichletBC(Q, u0, OB)
        
    Cinit = Constant(core_comp)
    C0.interpolate(Cinit)

    
    # Create blank arrays to append to 
    times = np.empty([0])
    core_mg = np.empty([0])

    # Timestepping
    i = 0
    zz = 0   
    while i < t:
    
        u_eval = C0(*mm) if mm_distance < DOLFIN_EPS else None
        computed_u = mesh.mpi_comm().gather(u_eval, root=0)
        computed_u = mesh.mpi_comm().bcast(computed_u, root=0)
        value = [u for u in computed_u if u is not None][0]
        times = np.append(times, i/(86400.0*365.0))
        core_mg = np.append(core_mg, value)
        
        solve(a==L, C0, Cbcs, solver_parameters={'linear_solver': 'cg', 'preconditioner': 'hypre_amg'}) #) #, s) solver_parameters={'linear_solver': 'mumps'}
        
        if mesh.mpi_comm().rank == 0:
            print(zz)
        
        i += dt 
        zz += 1

    return times, core_mg
      
# Intensive parameters used in modelling 

# Import files and meshes here

f_mesh = sys.argv[1]
f_dat = sys.argv[2]
f_dist = float(sys.argv[3])
f_temp = float(sys.argv[4])
element = sys.argv[5]
nts = float(sys.argv[6])
f_err = float(sys.argv[7])

#f_outfile = sys.argv[5]


##################################################
# Import mesh file from h5 file
mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), './{0}.h5'.format(f_mesh), "r")
hdf.read(mesh, "/mesh", False)
domains = MeshFunction('size_t', mesh, mesh.topology().dim())
hdf.read(domains, "/domains")
#boundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
#hdf.read(boundaries, "/boundaries")


# Define central point
mm = Point(0.0, 0.0, 0.0)
mm_cell, mm_distance = mesh.bounding_box_tree().compute_closest_entity(mm)


# Create file for input parameters? 

df_c = pd.read_csv(f_dat)

#Temp = df_c['T_C'].values[0] #1125
aSiO2 = df_c['aSiO2'].values[0]#0.71
Xan =  df_c['XAn'].values[0]# 0.87
A_x = df_c['A_{0}'.format(element)].values[0]

rim = df_c['rim_{0}'.format(element)].values[0] 
core = df_c['core_{0}'.format(element)].values[0] 

# Calculate initial plagioclase Mg content

# total time needs to be determined by crystal size and temperature


# Determine diffusion coefficient values based on element
Tx = f_temp + 273.0
Tx_i = 1.0/Tx

if element == 'Mg':
    Dtest = np.array([-5.45, -7.983, -3.54e4]) # Van Orman
if element == 'Ba': 
    Dtest = np.array([-12.32, -3.287, -3.995e4])
if element == 'Sr':
    Dtest = np.array([-12.81, -5.712, -3.239e4])


t_secs = ((f_dist**2.0)/D1(Tx_i, Xan, Dtest[0], Dtest[1], Dtest[2]))/3.0

t = t_secs/(86400.0*365.0)

#print(t)


equi_t = np.empty([0])


c_comp_t = mod_diff(t, f_temp, Xan, A_x, rim, core, element) 

tt = c_comp_t[0]
core_mg = c_comp_t[1]

count = 0
switch = False
for z in core_mg:
    if switch == False:
        if z < rim + f_err:
             #print(z)
             #print(count)
             #print(tt[count])
             equi_t = np.append(equi_t, tt[count])
             switch = True
    count += 1


# Input final dataframe

df = pd.read_csv('./Plag_3D_equi_time.csv')

df1 = pd.DataFrame({'Mesh':f_mesh, 'c half axis (um)':float(f_dist), 'Time (years)':tt,'Core_{0} (ppm)'.format(element):core_mg})
# Output to dataframe


# Append new rows to dataframe - should have identical headings - use append or concat?
df2 = pd.DataFrame({'Mesh':f_mesh, 'c half axis (um)':float(f_dist), 'Temperature (C)': float(f_temp), '{0} Equilibration Time (years)'.format(element):equi_t})

df3 = pd.concat([df, df2])

df1.to_csv('{0}_{1}_T{2}_3D_plag_model_output.csv'.format(f_mesh, f_dist, f_temp), sep=',')
df3.to_csv('./Plag_3D_equi_time.csv', sep=',')

#==============================================================================
# print(mesh.topology().dim())
# print(mesh.coordinates())
# print(type(mesh.coordinates()))
# 
# xdmf.write(mesh, XDMFFile.Encoding_HDF5)
#==============================================================================

#/bin/bash

# List of plagioclase z axis lengths used to scale mesh.

 

# output file names 

Pl_z="50.0"
f_temp="675.0"

file_name="high_albite_mesh"
f_dat=Krakatau_conditions.csv
element="Mg"
nts=1000.0
res=50.0
err=10.0

for x in $Pl_z ; do

echo ${x}

python3 albite_mesh_gen.py $x $file_name $res

# Use gmsh to create mesh

gmsh -3 "${file_name}.geo"

# convert mesh for dolfin 

dolfin-convert "${file_name}.msh" "${file_name}.xml"

python3 write_h5_mesh.py ${file_name}


for T in $f_temp ; do

echo ${T}

# Test mesh

mpirun -np 17 python3 3D_plag_diffusion_complex_2023_2.py ${file_name} ${f_dat} ${x} ${T} ${element} ${nts} ${err}

done

done

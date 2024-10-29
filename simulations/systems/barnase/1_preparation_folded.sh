#!/bin/bash

if [ -z "$ROOT" ]; then
    echo "ROOT variable is not defined. Exiting."
    exit 1
fi

if [ -z "$NAME" ]; then
    echo "NAME variable is not defined. Exiting."
    exit 1
fi

seed=$$

folder="$ROOT/systems/$NAME"
mkdir -p $folder

################################################################################
# MINIMIZATION IN VACUUM

folder_top="$folder/1-creating-topology"
mkdir -p "$folder_top"
pushd "$folder_top"

cp "$FF" .

## Generating the topology and the structure files, and setting up the force-field, 
# the water model and the histidine protonation state for the simulation:
{ yes 1; } | \
$gmx pdb2gmx -f "$PDB" -ff charmm36-jul2022 -water tip3p -ignh -his

popd

################################################################################
#CREATING THE WATER BOX AND SOLVATING THE PROTEIN

folder_water="$folder/2-adding-water"
mkdir -p $folder_water
pushd "$folder_water"

cp "$folder_top"/topol.top .
cp "$folder_top"/*.itp .
mv "$folder_top"/*.ff .

# Preparing the system by generating a box
$gmx editconf -f "$folder_top/conf.gro" -o newbox.gro -box 14 -bt dodecahedron \
     -c -princ > editconf.out <<EOF
1
EOF

# Add water molecules to the system
$gmx solvate -cp newbox.gro -cs spc216.gro -o water.gro -p topol.top

popd

################################################################################
#ADDING COUNTERIONS TO NEUTRALIZE THE SYSTEM

folder_ions="$folder/3-adding-ions"
mkdir -p $folder_ions

pushd "$folder_ions"

cp "$folder_water"/topol.top .
cp "$folder_water"/*.itp .
mv "$folder_water"/*.ff .

# Create a new TPR file for the system solvated
$gmx grompp -f "$ROOT/mdps/vac-minim.mdp" -c "$folder_water/water.gro" \
            -p topol.top -o solvated.tpr -maxwarn 10

# Neutralization step. Ions are added to the system and water molecules are 
# substituted by
$gmx genion -s solvated.tpr -o neut.gro -p topol.top -pname NA \
            -nname CL -conc 0.003 -neutral -seed $seed <<EOF
13
EOF

# Convert the generated gro into pdb
$gmx trjconv -f  neut.gro -s neut.gro -o neut.pdb -pbc mol -ur compact -center <<EOF
1
0
EOF

popd

########################################################################################################
#MINIMIZATION OF THE SOLVATED AND NEUTRALIZED SYSTEM

folder_minimization="$folder/4-solvated-minimization"
mkdir -p $folder_minimization

pushd "$folder_minimization"

# Minimization step
$gmx grompp -f "$folder/mdps/pbc-solv-minim.mdp" \
            -c "$folder_ions/neut.gro" -p "$folder_ions/topol.top" \
            -o "EM-neut.tpr" -maxwarn 10 >> EM-tpr.log

$gmx mdrun -deffnm EM-neut -nice 19

$gmx energy -f EM-neut.edr -o EM-neut_potential-energy.xvg <<EOF
12
0
EOF

$gmx trjconv -f EM-neut.gro -s EM-neut.tpr -o EM-neut.pdb \
             -pbc mol -ur compact -center <<EOF
1
0
EOF

popd

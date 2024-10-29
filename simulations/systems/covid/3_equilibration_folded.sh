#!/bin/bash

if [ -z "$ROOT" ]; then
    echo "ROOT variable is not defined. Exiting."
    exit 1
fi

if [ -z "$NAME" ]; then
    echo "NAME variable is not defined. Exiting."
    exit 1
fi

temp=298

################## MAIN BODY #######################

folder="$ROOT/systems/$NAME"
folder_ions="$folder/3-adding-ions"
folder_heating="$folder/5-heating"
folder_equil="$folder/6-equilibration"

mkdir -p "$folder_equil"
pushd "$folder_equil"

mdp1=md_equil1.mdp
mdp2=md_equil2.mdp
mdp3=md_equil3.mdp

cp "$folder/mdps/$mdp1" .
cp "$folder/mdps/$mdp2" .
cp "$folder/mdps/$mdp3" .

sed -i "s/^ref_t.*/ref_t = $temp/" "$mdp1"
sed -i "s/^ref_t.*/ref_t = $temp $temp/" "$mdp2"
sed -i "s/^ref_t.*/ref_t = $temp $temp/" "$mdp3"
sed -i "s/^gen_temp.*/gen_temp = $temp/" "$mdp1"
sed -i "s/^gen_temp.*/gen_temp = $temp/" "$mdp2"
sed -i "s/^gen_temp.*/gen_temp = $temp/" "$mdp3"

echo -n "The simulation is running in: " > RUNNING_INFORMATION.info
echo $HOSTNAME >> RUNNING_INFORMATION.info
echo "Initial time: " >> RUNNING_INFORMATION.info
date >> RUNNING_INFORMATION.info

# Start the initial steps for equilibration of the system

## First equilibration step
$gmx grompp -f "$mdp1" -po md_equil1_out.mdp -c "$folder_heating/heating.gro" \
            -r "$folder_heating/heating.gro" -p "$folder_ions/topol.top" \
            -o equil1.tpr -maxwarn 10

$gmx mdrun -s equil1.tpr -x equil1.xtc -c equil1.gro -e equil1.edr -nice 19 -dlb yes

## Second equilibration step
$gmx grompp -f "$mdp2" -po md_equil2_out.mdp -c equil1.gro -r equil1.gro \
            -t state.cpt -p "$folder_ions/topol.top" -o equil2.tpr -maxwarn 10

$gmx mdrun -s equil2.tpr -x equil2.xtc -c equil2.gro -e equil2.edr -nice 19 -dlb yes

## Third equilibration step
$gmx grompp -f "$mdp3" -po md_equil3_out.mdp -c equil2.gro -r equil2.gro \
            -t state.cpt -p "$folder_ions/topol.top" -o equil3.tpr -maxwarn 10

$gmx mdrun -s equil3.tpr -x equil3.xtc -c equil3.gro -e equil3.edr -nice 19 -dlb yes

$gmx trjconv -f equil3.xtc -s equil3.tpr -o equil3_preprocessed.xtc -pbc mol -ur compact -center <<EOF
1
0
EOF

$gmx trjconv -f equil3_preprocessed.xtc -s equil3.tpr -o equil3_processed.xtc -fit rot+trans <<EOF
1
0
EOF

echo "Final time: " >> RUNNING_INFORMATION.info
date >> RUNNING_INFORMATION.info

popd

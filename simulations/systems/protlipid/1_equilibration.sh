#!/bin/bash

if [ -z "$ROOT" ]; then
    echo "ROOT variable is not defined. Exiting."
    exit 1
fi

if [ -z "$NAME" ]; then
    echo "NAME variable is not defined. Exiting."
    exit 1
fi

## Main variables
equi_prefix=step6.%d_equilibration

#replica=r${SLURM_ARRAY_TASK_ID}
folder="$ROOT/systems/$NAME"
folder_prep="$folder/1-preparation"
folder_equil="$folder/2-equilibration"

mkdir -p "$folder_equil"
pushd "$folder_equil"

cp "$folder/mdps/step6.0_equilibration.mdp" .
cp "$folder/mdps/step6.1_equilibration.mdp" .
cp "$folder/mdps/step6.2_equilibration.mdp" .
cp "$folder/mdps/step6.3_equilibration.mdp" .
cp "$folder/mdps/step6.4_equilibration.mdp" .
cp "$folder/mdps/step6.5_equilibration.mdp" .
cp "$folder/mdps/step6.6_equilibration.mdp" .

cp "$folder_prep"/step6.0_minimization.gro .
cp "$folder_prep"/step6.0_minimization.tpr .

echo -n "The simulation is running in: " > RUNNING_INFORMATION.info
echo $HOSTNAME >> RUNNING_INFORMATION.info
echo "Initial time: " >> RUNNING_INFORMATION.info
date >> RUNNING_INFORMATION.info

# Equilibration
for cnt in $(seq 1 1 6); do
    pcnt=$((cnt - 1))
    istep=$(printf ${equi_prefix} ${cnt})
    pstep=$(printf ${equi_prefix} ${pcnt})
    if [ $cnt -eq 1 ]; then
        pstep="step6.0_minimization"
    fi

    srun $GMX_MPI_MOD grompp -f ${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ${pstep}.gro -p "$folder_prep/topol.top" -n "$folder_prep/index.ndx"
    srun $GMX_MPI_MOD mdrun -ntomp $OMP_NUM_THREADS -deffnm ${istep}
done

## Trajectory post-processing
srun $GMX_MPI_MOD trjconv -f ${istep}.xtc -s ${istep}.tpr -o ${istep}_preprocessed.xtc -pbc nojump <<EOF
0
EOF

srun $GMX_MPI_MOD trjconv -f ${istep}_preprocessed.xtc -s ${istep}.tpr -o ${istep}_processed.xtc -pbc mol -ur compact -center <<EOF
1
0
EOF

echo "Final time: " >> RUNNING_INFORMATION.info
date >> RUNNING_INFORMATION.info

popd



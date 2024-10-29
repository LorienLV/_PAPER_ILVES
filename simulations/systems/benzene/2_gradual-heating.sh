#!/bin/sh

if [ -z "$ROOT" ]; then
        echo "ROOT variable is not defined. Exiting."
        exit 1
fi

if [ -z "$NAME" ]; then
        echo "NAME variable is not defined. Exiting."
        exit 1
fi

## Main variables
temp=310

#replica=r${SLURM_ARRAY_TASK_ID}

## Heating (NVT step) ####
folder="$ROOT/systems/$NAME"
folder_prep="$folder/1-preparation"
folder_heating="$folder/2-heating"

mkdir -p "$folder_heating"
pushd "$folder_heating"

mdp="md_heat.mdp"
cp "$folder/mdps/$mdp" .

seed=$$
init_t=$(($temp - 300))

sed -i "s/^gen_seed.*/gen_seed =  $seed/" "$mdp"
sed -i "s/^ref_t.*/ref_t = $init_t/" "$mdp"
sed -i "s/^; Generate velocites.*/; Generate velocites is on at $init_t K/" "$mdp"
sed -i "s/^gen_temp.*/gen_temp =  $init_t/" "$mdp"
sed -i "s/^gen_seed.*/gen_seed =  $seed/" "$mdp"

echo -n "The simulation is running in: " >RUNNING_INFORMATION.info
echo $HOSTNAME >>RUNNING_INFORMATION.info
echo "Initial time: " >>RUNNING_INFORMATION.info
date >>RUNNING_INFORMATION.info

# Loop for increasing the the system's T
srun $GMX_MPI_MOD grompp -f "$mdp" -po md_heat-0-out.mdp -c "$folder_prep/$NAME-min.gro" \
        -r "$folder_prep/$NAME-min.gro" -p "$folder_prep/$NAME.top" \
        -o heating-0.tpr -maxwarn 10

srun $GMX_MPI_MOD mdrun -s heating-0.tpr -c heating.gro -e heating.edr -x heating.xtc \
        -g heating.log -nice 19

for i in {1..6}; do
        new_t=$(($init_t + $i * 50))
        new_ext=$((25000 + $i * 25000))

        sed -i "s/^ref_t.*/ref_t = $new_t/" "$mdp"
        sed -i "s/^; Generate velocites.*/; Generate velocites is off at $new_t K/" "$mdp"
        sed -i "s/^gen_vel.*/gen_vel = no/" "$mdp"
        sed -i "s/^gen_temp.*/gen_temp = $new_t/" "$mdp"
        sed -i "s/^nsteps.*/nsteps = $new_ext/" "$mdp"
        sed -i "s/^gen_seed.*/gen_seed = $seed/" "$mdp"

        srun $GMX_MPI_MOD grompp -f "$mdp" -po md_heat-$i-out.mdp -c heating.gro -r heating.gro \
                -p "$folder_prep/$NAME.top" -o heating-$i.tpr -maxwarn 10

        srun $GMX_MPI_MOD mdrun -s heating-$i.tpr -c heating.gro -cpi state.cpt -e heating.edr \
                -x heating.xtc -g heating.log -nice 19
done

echo "Final time: " >>RUNNING_INFORMATION.info
date >>RUNNING_INFORMATION.info

popd

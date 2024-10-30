#!/bin/bash

# Perform production runs using different constraint solvers, tolerances and configurations.
# $1 -> output folder
# $2 -> all-bonds, h-bonds or h-angles

# NTASKS, NTASKS_NODE, NTASKS_SOCKET, NTHREADS_PER_TASKS,
configs=("1,1,1,1" "1,1,1,56")

# If you are using single precision, the tolerance must be lower than 0.000001.
tols=(0.0001 0.00000001 0.000000000001)

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters"
    echo "Usage: prod_run.sh <output_folder> <all-bonds|h-bonds|h-angles>"
    exit 1
fi

root_folder="$1"
if [ -z "$root_folder" ]; then
    echo "Output folder is not defined. Exiting."
    exit 1
fi

case "$2" in
    all-bonds)
        constraints="all-bonds"
        ;;
    h-bonds)
        constraints="h-bonds"
        ;;
    h-angles)
        constraints="h-angles"
        ;;
    *)
        echo "Invalid constraint type. Exiting."
        exit 1
        ;;
esac

if [ -z "$SIMULATIONS_DIR" ]; then
    echo "ROOT variable is not defined. Exiting."
    exit 1
fi

if [ -z "$NAME" ]; then
    echo "NAME variable is not defined. Exiting."
    exit 1
fi

grofile="$SIMULATIONS_DIR/systems/$NAME/files-for-prod/gro.gro"
cptfile="$SIMULATIONS_DIR/systems/$NAME/files-for-prod/cpt.cpt"
topolfile="$SIMULATIONS_DIR/systems/$NAME/files-for-prod/top.top"
ndxfile="$SIMULATIONS_DIR/systems/$NAME/files-for-prod/ndx.ndx"

#
# REMEMBER TO SET:
# gen_vel      = no
# gen-seed     = X
# ld-seed      = X
# TO ENSURE REPRODUCIBILITY
#
mdmdp="$SIMULATIONS_DIR/systems/$NAME/mdps/md_prod.mdp"

export GMX_MAXBACKUP=-1
export CONSTRAINTS_VERBOSE=0

export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Run GROMACS using a constraint solver.
# $1 = gmx binary
# $2 = solver folder
# $3 = tpr file
# $4 = config
run_solver() {
    gmx_bin="$1"
    solver_folder="$2"
    tpr="$3"
    config="$4"

    ntasks=$(echo $config | cut -d',' -f1)
    ntasks_node=$(echo $config | cut -d',' -f2)
    ntasks_socket=$(echo $config | cut -d',' -f3)
    nthreads_task=$(echo $config | cut -d',' -f4)

    echo "        $solver_folder..."

    mkdir "$solver_folder" &>/dev/null
    pushd "$solver_folder" &>/dev/null

    sched_file="run.sh"

    echo "#!/bin/bash" > "$sched_file"

    # echo "" >> "$sched_file"
    # echo "#SBATCH --ntasks=$ntasks" >> "$sched_file"
    # echo "#SBATCH --ntasks-per-node=$ntasks_node" >> "$sched_file"
    # echo "#SBATCH --ntasks-per-socket=$ntasks_socket" >> "$sched_file"
    # echo "#SBATCH --cpus-per-task=$nthreads_task" >> "$sched_file"
    # echo "#SBATCH --job-name=$root_folder-$config" >> "$sched_file"

    echo "" >> "$sched_file"
    # echo "export SRUN_CPUS_PER_TASK=$nthreads_task" >> "$sched_file"
    echo "export OMP_NUM_THREADS=$nthreads_task" >> "$sched_file"
    echo "export OMP_PROC_BIND=true" >> "$sched_file"
    echo "export OMP_PLACES=cores" >> "$sched_file"

    command="$gmx_bin mdrun -s "../../$tpr" -noconfout \
             -g gromacs.log -ntomp $nthreads_task 1>out.out 2>err.err"

    echo "" >> "$sched_file"
    echo "$command" >> "$sched_file"

    bash "$sched_file"

    popd &>/dev/null
}

# Generate mdp file based on an old mdp file and a string of parameters.
# $1 = old mdp file
# $2 = new mdp file
# $3 = new parameters
generate_mdp() {
    old_mdp="$1"
    new_mdp="$2"
    new_parameters="$3"

    cp "$old_mdp" "$new_mdp"

    echo "; Added by the script" >> "$new_mdp"

    new_parameters+=""$'\n'"constraints=$constraints"

    for line in $new_parameters; do
        trilled_line=$(echo "$line" | sed 's/\s//g')
        option=$(echo "$trilled_line" | cut -d'=' -f1)
        value=$(echo "$trilled_line" | cut -d'=' -f2)

        # The option can use both - and _ as separators.
        option_regex=$(echo "$option" | sed -r 's/[-_]/\[-_\]/g')

        # Comment the option if it already exists.
        sed -i -E "s/^[[:space:]]*${option_regex}[[:space:]]*=.*$/;&/gI" "$new_mdp"
        # Append the new option at the end of the new mdp file.
        echo "${option}=${value}" >> "$new_mdp"
    done
}

# Generate trp file based on an mdp file and a string of parameters.
# $1 = gmx binary
# $2 = old mdp file
# $3 = new mdp file
# $4 = new tpr file
# $5 = new parameters
generate_tpr() {
    gmx_bin="$1"
    old_mdp="$2"
    new_mdp="$3"
    new_tpr="$4"
    new_parameters="$5"

    echo "    $new_tpr..."

    generate_mdp "$old_mdp" "$new_mdp" "$new_parameters"

    command="$gmx_bin grompp -f "$new_mdp""
    command+=" -c "$grofile" -r "$grofile" -t "$cptfile" -p "$topolfile""

    if [ -f "$ndxfile" ]; then
        command+=" -n "$ndxfile""
    fi

    command+=" -o "$new_tpr""
    command+=" -maxwarn 2"

    # Create the tpr file
    eval "$command" &> "${new_tpr}.txt"
}

for tol in ${tols[@]}; do

    echo "Tolerance: $tol"

    tol_folder="${root_folder}/t${tol}"
    mkdir -p "$tol_folder" &>/dev/null
    pushd "$tol_folder" &>/dev/null

    echo "Generating tpr files..."

    shake_mdp="shake.mdp"
    plincs_mdp="plincs_i2o4.mdp"
    mplincs_o4_mdp="mplincs_o4.mdp"
    mplincs_o8_mdp="mplincs_o8.mdp"
    ilves_mdp="ilves.mdp"
    ilvesf_mdp="ilvesf.mdp"

    shake_tpr="shake.tpr"
    plincs_tpr="plincs_i2o4.tpr"
    mplincs_o4_tpr="mplincs_o4.tpr"
    mplincs_o8_tpr="mplincs_o8.tpr"
    ilves_tpr="ilves.tpr"
    ilvesf_tpr="ilvesf.tpr"

    # SHAKE
    if true; then
        generate_tpr "$GMX_MPI_MOD" "$mdmdp" "$shake_mdp" "$shake_tpr" \
                     "constraint-algorithm=shake
                      shake-tol=${tol}"
    fi

    # PLINCS
    if true; then
        generate_tpr "$GMX" "$mdmdp" "$plincs_mdp" "$plincs_tpr" \
                     "constraint-algorithm=lincs
                      lincs-iter=2
                      lincs-order=4"
    fi

    # MPLINCS-O4
    if true; then
        generate_tpr "$GMX_MPI_MOD" "$mdmdp" "$mplincs_o4_mdp" "$mplincs_o4_tpr" \
                     "constraint-algorithm=lincs
                      lincs-tol=${tol}
                      lincs-iter=1000
                      lincs-order=4"
    fi

    # MPLINCS-O8
    if true; then
        generate_tpr "$GMX_MPI_MOD" "$mdmdp" "$mplincs_o8_mdp" "$mplincs_o8_tpr" \
                     "constraint-algorithm=lincs
                      lincs-tol=${tol}
                      lincs-iter=1000
                      lincs-order=8"
    fi

    # ILVES
    if true; then
        generate_tpr "$GMX_MPI_MOD" "$mdmdp" "$ilves_mdp" "$ilves_tpr" \
                     "constraint-algorithm=ilves
                      ilves-tol=${tol}
                      ilves-mpi-neigh=4"
    fi

    # ILVESF
    if true; then
        generate_tpr "$GMX_MPI_MOD" "$mdmdp" "$ilvesf_mdp" "$ilvesf_tpr" \
                     "constraint-algorithm=ilvesf
                      ilves-tol=${tol}
                      ilves-mpi-neigh=4"
    fi

    echo "Running simulations..."

    for config in ${configs[@]}; do
        ntasks=$(echo $config | cut -d',' -f1)
        ntasks_node=$(echo $config | cut -d',' -f2)
        ntasks_socket=$(echo $config | cut -d',' -f3)
        nthreads_task=$(echo $config | cut -d',' -f4)

        config_folder="t${ntasks}n${ntasks_node}s${ntasks_socket}c${nthreads_task}"
        mkdir -p "$config_folder" &>/dev/null
        pushd "$config_folder" &>/dev/null

        echo "    Config:"
        echo "        Tasks: ${ntasks}"
        echo "        Tasks per Node: ${ntasks_node}"
        echo "        Tasks per Socket: ${ntasks_socket}"
        echo "        Cores per Task: ${nthreads_task}"

        echo "    Simulation:"

        #
        # SHAKE
        #
        if true; then
            if [[ $ntasks -eq 1 ]]; then
                run_solver "$GMX_MPI_MOD" "SHAKE" "$shake_tpr" "$config"
            fi
        fi

        #
        # PLINCS
        #
        if true; then
            run_solver "$GMX" "P-LINCS-I2O4" "$lincs_tpr" "$config"
        fi

        #
        # MPLINCS-O4
        #
        if true; then
            run_solver "$GMX_MPI_MOD" "P-LINCS-O4" "$lincs_o4_tpr" "$config"
        fi

        #
        # MPLINCS-O8
        #
        if true; then
            run_solver "$GMX_MPI_MOD" "P-LINCS-O8" "$lincs_o8_tpr" "$config"
        fi

        #
        # ILVES
        #
        if true; then
            run_solver "$GMX_MPI_MOD" "ILVES-ASYM" "$ilves_asym_tpr" "$config"
        fi

        #
        # ILVESF
        #
        if true; then
            run_solver "$GMX_MPI_MOD" "ILVES-SYM" "$ilves_sym_tpr" "$config"
        fi

        popd &>/dev/null # run folder
    done
    popd &>/dev/null # tpr folder
done
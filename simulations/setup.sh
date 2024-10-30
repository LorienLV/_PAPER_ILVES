#!/bin/bash

# The parent folder of this script.
export SIMULATIONS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export GMX_MPI="PATH_TO_GROMACS_DOUBLE_INSTALLATION/gmx_mpi_d"
export GMX_MPI_MOD="PATH_TO_GROMACS_WITH_ILVES_DOUBLE_INSTALLATION/gmx_mpi_d"

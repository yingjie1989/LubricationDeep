#!/bin/bash
#
#SBATCH --job-name=praft
#SBATCH --output=praft_output.txt
#SBATCH --error=praft_error.txt
#SBATCH --ntasks=32
#SBATCH --time=24:00:00





mpiexec ../../Moose_work_praft_model2_newAC/praft_fracture-opt -i 2D_TOTAL.i

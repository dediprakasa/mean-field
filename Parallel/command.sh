#!/bin/bash

#SBATCH -N 8
#SBATCH -n 32
#SBATCH -p hp32c

path='/3DMF'
prefix='3DMF'

echo "Program starts at:"
date '+%A %D %X'

echo "Running test"

mpif90 -g -o 3DMF.exe global.f90 allocation.f90 deallocation.f90 main.f90 fermi.f90 rtbis.f90 Optical_Response.f90 -llapack

mpirun -np $SLURM_NTASKS ./3DMF.o < input_main > logfile.dat

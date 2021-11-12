#!/bin/bash
 
#SBATCH -n 1                       # number of cores
#SBATCH -t 0-2:00                   # wall time (D-HH:MM)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)

#module purge
#module load matlab/2020a
#module load gcc/9.2.0
module load intel/2020.2

#make -f Makefile.intel
./intel_HQRRP_single 496 10201 mdp30r496c10201.bin pf30.bin  

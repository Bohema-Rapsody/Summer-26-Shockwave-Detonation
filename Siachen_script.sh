#!/bin/bash
#SBATCH --job-name=Cu_30nm_equilib
#SBATCH --time=48:00:00
#SBATCH --mem=26000
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=av704@cam.ac.uk

cd ~/Documents/Summer_Research_2026/LAMMPS/Summer-26-Shockwave-Detonation
unset DISPLAY

export OMP_NUM_THREADS=1

mpirun -np 24 ~/lammps/build/lmp -in Cu-N2/Cu-N2.in
#!/bin/bash
#SBATCH -p cclake
#SBATCH -A KATERIS-SL3-CPU
#SBATCH --job-name=N2_shock_sim
#SBATCH --time=03:00:00
#SBATCH --mem=3410
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=56
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=av704@cam.ac.uk

cd ~/Documents/Summer_Research_2026/LAMMPS/Summer-26-Shockwave-Detonation
unset DISPLAY

export OMP_NUM_THREADS=1

srun ~/lammps/build/lmp -in Cu-N2/N2_Cu_combine.in
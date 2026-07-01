# Summer-26-Shockwave-Detonation
LAMMPS files for simulation


Running on Siachen:

#!/bin/bash
#SBATCH --job-name=test_Argon
#SBATCH --time=00:05:00
#SBATCH --mem=26000
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=av704@cam.ac.uk

cd ~/Documents/Summer_Research_2026/LAMMPS/Example/Ar
unset DISPLAY

export OMP_NUM_THREADS=1

mpirun -np 24 ~/lammps/build/lmp -in ar.in



initial modelling task - getting used to LAMMPS:
model the PVT relations of Ar and N2 gas at different temperatures, pressures and densities to determine at which point ideal laws fail
Both gases uses a potential field model using Lennard-Jones - pairwise calculation with some cut-off radius
 - see parameter properties epsilon and sigma through empirical fitting research

metal structure - copper:
creating embedded-atom model (EAM) for metal lattice potential fields
Cu Potentials: https://www.ctcms.nist.gov/potentials/system/Cu/

Chosen potential: https://www.ctcms.nist.gov/potentials/entry/2018--Etesami-S-A-Asadi-E--Cu/

using Lennard-jones potential fields to model the mixing of N2 and Cu particles

measurements and comparisons:
heating from shockwave compared to thermal conduction
comparing solid particle heating to liquid heating - comparing viscous forces and dissipation contribution to heating
comparing drag measurements to standards Stoke's drag models - momentum and energy transfer are ultimately coupled
Matching non-dimensional values to achieve similar results to expected situations in practise
 - Knudsen numbers need to match as the real droplets will be large compared to the mean free path of air molecules


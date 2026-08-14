# Summer-26-Shockwave-Detonation
LAMMPS files for simulation


Running on Siachen:

#!/bin/bash
#SBATCH --job-name=N2_shock_sim
#SBATCH --time=00:10:00
#SBATCH --mem=26000
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=av704@cam.ac.uk

cd ~/Documents/Summer_Research_2026/LAMMPS/Summer-26-Shockwave-Detonation
unset DISPLAY

export OMP_NUM_THREADS=1

mpirun -np 24 ~/lammps/build/lmp -in N2_shock/N2_shock_contd.in


Installing LAMMPS on CSD3:

module purge
module load rhel8/default-icl
 
git clone https://github.com/lammps/lammps.git
cd lammps
 
mkdir build
cd build
 
cmake ../cmake \
    -D BUILD_MPI=on \
    -D PKG_MOLECULE=on \
    -D PKG_MEAM=on \
    -D PKG_MANYBODY=on
 
make -j4

check for packages:
mpiexec -n 5 ./lmp -h


Running on CSD3 (av704@login.hpc.cam.ac.uk):

#!/bin/bash
#SBATCH -p cclake
#SBATCH -A KATERIS-SL3-CPU
#SBATCH --job-name=N2_shock_sim
#SBATCH --time=03:00:00
#SBATCH --mem=26000
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=56
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=av704@cam.ac.uk

cd ~/Documents/Summer_Research_2026/LAMMPS/Summer-26-Shockwave-Detonation
unset DISPLAY

export OMP_NUM_THREADS=1

srun ~/lammps/build/lmp -in Cu-N2/N2_Cu_combine.in


...

Useful slurm commands:

squeue -u av704
scontrol show job ########


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



SSH commands:

Copying files:

To Gate (from PC project directory):
scp -r ./Summer-26-Shockwave-Detonation av704@gate.eng.cam.ac.uk:"~/Documents/Summer_Research_2026/LAMMPS"

scp -r ./* av704@gate.eng.cam.ac.uk:"~/Documents/Summer_Research_2026/LAMMPS/Summer-26-Shockwave-Detonation"



To Siachen (from Gate):
rsync -ruvP ~/Documents/Summer_Research_2026/LAMMPS/Summer-26-Shockwave-Detonation av704@siachen.eng.cam.ac.uk:"~/Documents/Summer_Research_2026/LAMMPS"

rsync -ruvP ~/Documents/Summer_Research_2026/LAMMPS/Summer-26-Shockwave-Detonation/* av704@siachen.eng.cam.ac.uk:"~/Documents/Summer_Research_2026/LAMMPS/Summer-26-Shockwave-Detonation"



To Gate (from Siachen via Gate):
rsync -ruvP av704@siachen.eng.cam.ac.uk:"~/Documents/Summer_Research_2026/LAMMPS/Summer-26-Shockwave-Detonation" ~/Documents/Summer_Research_2026/LAMMPS

To PC (from Gate via PC project directory):
scp -r av704@gate.eng.cam.ac.uk:"~/Documents/Summer_Research_2026/LAMMPS/Summer-26-Shockwave-Detonation/*" ./

To Laptop (from PC via Laptop project directory)
...

To CSD3 (from PC project directory):
scp -r ./Summer-26-Shockwave-Detonation av704@login.hpc.cam.ac.uk:"~/Documents/Summer_Research_2026/LAMMPS"
scp -r ./* av704@login.hpc.cam.ac.uk:"~/Documents/Summer_Research_2026/LAMMPS/icelake"

To PC (from CSD3 via PC project dir):
scp -r av704@login.hpc.cam.ac.uk:"~/Documents/Summer_Research_2026/LAMMPS" ./Summer-26-Shockwave-Detonation
scp -r av704@login.hpc.cam.ac.uk:"~/Documents/Summer_Research_2026/LAMMPS/icelake/*" ./
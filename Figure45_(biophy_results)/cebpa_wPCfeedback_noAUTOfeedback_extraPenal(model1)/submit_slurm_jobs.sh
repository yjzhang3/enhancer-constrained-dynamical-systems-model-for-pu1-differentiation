#!/bin/bash

#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=05:00:00
#SBATCH --mail-user=yz183@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --array=200-250

module purge
module load MATLAB/2024b
which matlab

srun matlab -nodisplay -r "ODE_PS_NOTS_Main(0,[2,2],'N'); exit;"

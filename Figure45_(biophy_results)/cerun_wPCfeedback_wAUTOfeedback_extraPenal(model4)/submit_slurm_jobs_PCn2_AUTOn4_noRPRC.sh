#!/bin/bash

#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=08:00:00
#SBATCH --mail-user=yz183@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --array=300-350

#1,3,4,8,9,10,11,12,16,17,18,19,20

module purge
module load MATLAB/2024b
which matlab

srun matlab -nodisplay -r "ODE_PS_NOTS_Main([2,2],[4,4],[0,0],'Y',0); exit;"


#!/bin/bash

#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=03:00:00
#SBATCH --mail-user=yz183@rice.edu
#SBATCH --mail-type=ALL
#SBATCH --array=200-250

module purge
module load MATLAB/2024b
which matlab

#srun matlab -nodisplay -r "debug_model2_v3(0); exit;"
srun matlab -nodisplay -r "debug_model2_v4(2,0,'model2'); exit;"

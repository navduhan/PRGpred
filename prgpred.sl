#!/bin/bash

#SBATCH -t 1-0:00       # hours:minutes runlimit after which job will be killed
#SBATCH --nodelist chela-g01
#SBATCH --partition mahaguru
#SBATCH --job-name rslpred2       # Job name
#SBATCH -o /home/naveen/PRGpred/slurm/%j.out # File to which standard out will be written
#SBATCH -e /home/naveen/PRGpred/slurm/%j.err       # File to which standard err will be written
#SBATCH -W



module load ml-gpu
echo python /home/naveen/PRGpred/PRGpred.py -i $1 -l $2 -od $3
python /home/naveen/PRGpred/PRGpred.py -i $1 -l $2 -od $3

module purge

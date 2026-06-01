#!/bin/bash
#SBATCH --job-name=${INPUT}
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=gpu
#SBATCH --gres=gpu:p100:1
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
module load matlab
hostname 
matlab -nodisplay -batch "problemFileName = '${INPUT}'; main3; exit;"
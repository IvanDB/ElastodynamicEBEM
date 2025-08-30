#!/bin/sh
#SBATCH --job-name=launchP100std
#SBATCH --output=%x.output
#SBATCH --error=%x.errors
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=gpu
#SBATCH --gres=gpu:p100:1
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00

hostname
module load matlab/R2023b

matlab -nodisplay -batch "jobLauncher; exit;"
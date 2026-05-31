#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=128
#SBATCH --gres=gpu:p100:1


echo $COMMAND
matlab -nodisplay -batch $COMMAND

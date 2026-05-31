#!/bin/bash

module load cuda
module load matlab
matlab -nodisplay -batch "$COMMAND main3;"

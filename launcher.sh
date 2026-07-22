#!/bin/bash

module load cuda
module load matlab
matlab -nodisplay -batch "$COMMAND eebem.main;"

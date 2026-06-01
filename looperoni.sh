#!/bin/sh


inputPBs=("DesCop-cube" "barH1") #
meshTypes=("symm" "asym") #
levs=(1 2 3) # 

 
for inputPB in "${inputPBs[@]}"; do
for meshType in "${meshTypes[@]}"; do
for lev in "${levs[@]}"; do
 
 echo ${inputPB}-${meshType}_lev${lev}
 sbatch	--job-name="${inputPB}-${meshType}_lev${lev}" \
	--export=INPUT=input_${inputPB}-${meshType}_lev${lev}.txt \
	./main3.sh 2>/dev/null

done
done
done

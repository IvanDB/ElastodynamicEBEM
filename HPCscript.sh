#!/bin/bash

basePath=$(pwd)
source "$basePath/parameterLists.config"

baseTemplate="basePath=\"$basePath\"; inputStruct.problemName=[[PB]]; inputStruct.meshType=[[MT]]; inputStruct.meshLevel=[[ML]]; inputStruct.beta=[[beta]];"
quadTemplate="inputStruct.numSRext=[[Se]]; inputStruct.numGHext=[[Ge]]; inputStruct.numSRint=[[Si]]; inputStruct.numGHint==[[Gi]]; inputStruct.numSNGLR==[[SG]]; inputStruct.numBOUND=[[BD]];"

if [ ${#pbNames[@]} -eq 0 ]; then
    echo "Select at least one problem"
    exit 0
fi

for pbName in "${pbNames[@]}"; do
for meshType in "${meshTypes[@]:-"symm"}"; do
for meshLvl in "${meshLvls[@]:-1}"; do
for betaVal in "${betaVals[@]:-1}"; do
   
    baseData=$(sed  -e "s/\[\[PB\]\]/\"${pbName}\"/"    \
		            -e "s/\[\[MT\]\]/\"${meshType}\"/"  \
		            -e "s/\[\[ML\]\]/${meshLvl}/"       \
		            -e "s/\[\[beta\]\]/${betaVal}/"     \
		            <<< "$baseTemplate")

    if [ ${#quadID[@]} -ne 0 ]; then
        for quadID in "${quadID[@]}"; do
            cmdString="${baseData} inputStruct.quadID=$quadID;"
            cmdBuffer+=("$cmdString")

            jobName="${pbName}-${meshType}_${meshLvl}-${betaVal}_\'ID$quadID\'"
            nameBuffer+=("$jobName")
        done
        continue
    fi

    for numSRext in "${numsSRext[@]:-16}"; do
    for numGHext in "${numsSRext[@]:-1}"; do
    for numSRint in "${numsSRint[@]:-64}"; do
    for numGHint in "${numsSRint[@]:-3}"; do
    for numSNGLR in "${numsSNGLR[@]:-256}"; do
    for numBOUND in "${numsBOUND[@]:-256}"; do
        
        quadData=$(sed  -e "s/\[\[Se\]\]/${numSRext}/" \
                        -e "s/\[\[Ge\]\]/${numGHext}/" \
                        -e "s/\[\[Si\]\]/${numSRint}/" \
                        -e "s/\[\[Gi\]\]/${numGHint}/" \
                        -e "s/\[\[SG\]\]/${numSNGLR}/" \
                        -e "s/\[\[BD\]\]/${numBOUND}/" \
                        <<< "$quadTemplate")

        cmdString="${baseData} ${quadData}"
	    cmdBuffer+=("$cmdString")

        jobName="${pbName}-${meshType}_${meshLvl}-${betaVal}_'${numSRext}-${numGHext}_${numSRint}-${numGHint}_${numSNGLR}_${numBOUND}'"
        nameBuffer+=("$jobName")
    done
    done
    done
    done
    done
    done

done
done
done
done

if [ ${#cmdBuffer[@]} -eq 0 ]; then
    echo "Error in setup script: no command to execute."
    exit 1
fi

for i in "${!cmdBuffer[@]}"; do
    cmd="${cmdBuffer[i]}"
    name="${nameBuffer[i]}"

    echo "Submitting job $name"
    export COMMAND=${cmd}
    sbatch  --job-name=$name                \
            --output="logs/%j.log"          \
            --qos=gpu                       \     
            --partition=gpu                 \
            --nodes=1                       \
            --ntasks-per-node=$poolSize     \
            --mem=$RAMsize                  \
            --gres=gpu:$GPUtype:$GPUcount   \
            "./launcher.sh"
done
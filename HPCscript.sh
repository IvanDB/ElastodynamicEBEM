#!/bin/bash

basePath=$(pwd)
source "$basePath/parameterLists.config"

baseTemplate="basePath=\"$basePath\"; inputStruct.problemName=[[PB]]; [[MT]] [[ML]] [[beta]]"
quadTemplate="[[Se]] [[Ge]] [[Si]] [[Gi]] [[SG]] [[BD]]"

if [ ${#pbNames[@]} -eq 0 ]; then
    echo "Select at least one problem"
    exit 0
fi

for pbName in "${pbNames[@]}"; do
for meshType in "${meshTypes[@]:-""}"; do
for meshLvl in "${meshLvls[@]:--1}"; do
for betaVal in "${betaVals[@]:--1}"; do

    [[ $meshType == "" ]] && mt="" || mt="inputStruct.meshType=\"$meshType\";"
    [[ $meshLvl -eq -1 ]] && ml="" || ml="inputStruct.meshLevel=$meshLvl;"
    [[ $betaVal -eq -1 ]] && bv="" || bv="inputStruct.beta=$betaVal;"

    baseData=$(sed  -e "s/\[\[PB\]\]/\"${pbName}\"/"    \
		            -e "s/\[\[MT\]\]/$mt/"              \
		            -e "s/\[\[ML\]\]/${ml}/"            \
		            -e "s/\[\[beta\]\]/${bv}/"          \
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

    for numSRext in "${numsSRext[@]:--1}"; do
    for numGHext in "${numsGHext[@]:--1}"; do
    for numSRint in "${numsSRint[@]:--1}"; do
    for numGHint in "${numsGHint[@]:--1}"; do
    for numSNGLR in "${numsSNGLR[@]:--1}"; do
    for numBOUND in "${numsBOUND[@]:--1}"; do
        
        [[ $numSRext -eq -1 ]] && se="" || se="inputStruct.numSRext=$numSRext;"
        [[ $numGHext -eq -1 ]] && ge="" || ge="inputStruct.numGHext=$numGHext;"
        [[ $numSRint -eq -1 ]] && si="" || si="inputStruct.numSRint=$numSRint;"
        [[ $numGHint -eq -1 ]] && gi="" || gi="inputStruct.numGHint=$numGHint;"
        [[ $numSNGLR -eq -1 ]] && sg="" || sg="inputStruct.numSNGLR=$numSNGLR;"
        [[ $numBOUND -eq -1 ]] && bd="" || bd="inputStruct.numBOUND=$numBOUND;"

        quadData=$(sed  -e "s/\[\[Se\]\]/${se}/" \
                        -e "s/\[\[Ge\]\]/${ge}/" \
                        -e "s/\[\[Si\]\]/${si}/" \
                        -e "s/\[\[Gi\]\]/${gi}/" \
                        -e "s/\[\[SG\]\]/${sg}/" \
                        -e "s/\[\[BD\]\]/${bd}/" \
                        <<< "$quadTemplate")

        cmdString="${baseData} ${quadData}"
	    cmdBuffer+=("$cmdString")

        jobName="${pbName}-${mt}_${ml}-${bv}_'${se}-${ge}_${si}-${gi}_${sg}_${bd}'"
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
    echo $COMMAND
    # sbatch  --job-name=$name                \
    #         --output="logs/%j_%x.log"          \
    #         --qos=gpu                       \
    #         --partition=gpu                 \
    #         --nodes=1                       \
    #         --ntasks-per-node=$poolSize     \
    #         --mem=$RAMsize                  \
    #         --gres=gpu:$GPUtype:$GPUcount   \
    #         "./launcher.sh"
done
#!/bin/bash

basePath=$(pwd)
source "$basePath/parameterLists.config"

baseTemplate="basePath=\"$basePath\"; inputStruct.problemName=[[PB]]; [[MT]] [[ML]] [[BV]]"
quadTemplate="[[SE]] [[GE]] [[SI]] [[GI]] [[SG]] [[BD]]"

if [ ${#pbNames[@]} -eq 0 ]; then
    echo "Select at least one problem"
    exit 0
fi

for pbName in "${pbNames[@]}"; do
for meshType in "${meshTypes[@]:-""}"; do
for meshLvl in "${meshLvls[@]:--1}"; do
for betaVal in "${betaVals[@]:--1}"; do

    [[ $meshType == "" ]] && (mt=""; j_mt="") || (mt="inputStruct.meshType=\"$meshType\";"; j_mt="_$meshType")
    [[ $meshLvl -eq -1 ]] && (ml=""; j_ml="") || (ml="inputStruct.meshLevel=$meshLvl;"; j_ml="_$meshLvl")
    [[ $betaVal -eq -1 ]] && (bv=""; j_bv="") || (bv="inputStruct.beta=$betaVal;"; j_bv="_$betaVal")

    baseData=$(sed  -e "s/\[\[PB\]\]/\"${pbName}\"/"    \
		            -e "s/\[\[MT\]\]/$mt/"              \
		            -e "s/\[\[ML\]\]/$ml/"              \
		            -e "s/\[\[BV\]\]/$bv/"              \
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
        
        [[ $numSRext -eq -1 ]] && (se=""; j_se="") || (se="inputStruct.numSRext=$numSRext;"; j_se="_$numSRext")
        [[ $numGHext -eq -1 ]] && (ge=""; j_ge="") || (ge="inputStruct.numGHext=$numGHext;"; j_ge="_$numGHext")
        [[ $numSRint -eq -1 ]] && (si=""; j_si="") || (si="inputStruct.numSRint=$numSRint;"; j_si="_$numSRint")
        [[ $numGHint -eq -1 ]] && (gi=""; j_gi="") || (gi="inputStruct.numGHint=$numGHint;"; j_gi="_$numGHint")
        [[ $numSNGLR -eq -1 ]] && (sg=""; j_sg="") || (sg="inputStruct.numSNGLR=$numSNGLR;"; j_sg="_$numSNGLR")
        [[ $numBOUND -eq -1 ]] && (bd=""; j_bd="") || (bd="inputStruct.numBOUND=$numBOUND;"; j_bd="_$numBOUND")

        quadData=$(sed  -e "s/\[\[SE\]\]/$se/" \
                        -e "s/\[\[GE\]\]/$ge/" \
                        -e "s/\[\[SI\]\]/$si/" \
                        -e "s/\[\[GI\]\]/$gi/" \
                        -e "s/\[\[SG\]\]/$sg/" \
                        -e "s/\[\[BD\]\]/$bd/" \
                        <<< "$quadTemplate")

        cmdString="${baseData} ${quadData}"
	    cmdBuffer+=("$cmdString")

        jobName="${pbName}${j_mt}${j_ml}${j_bv}${j_se}${j_ge}${j_si}${j_gi}${j_sg}${j_bd}"
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
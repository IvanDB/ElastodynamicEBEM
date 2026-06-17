#!/bin/bash

basePath=$(pwd)
source "$basePath/parameterLists.config"

baseTemplate="basePath=\"$basePath\"; inputStruct.pbName=[[PB]]; inputStruct.formID=[[FI]]; [[MT]] [[ML]] [[BM]] [[TM]]"
quadTemplate="[[QT]] [[SE]] [[GE]] [[SI]] [[GI]] [[SG]] [[BD]]"

if [ ${#pbNames[@]} -eq 0 ]; then
    echo "Select at least one problem"
    exit 0
fi

if [ ${#formIDs[@]} -eq 0 ]; then
    echo "Select at least one formulation"
    exit 0
fi

for pbName in "${pbNames[@]}"; do
for formID in "${formIDs[@]}"; do
for meshType in "${meshTypes[@]:-""}"; do
for meshLvl in "${meshLvls[@]:--1}"; do
for betaMult in "${betaMults[@]:--1}"; do
for timeMult in "${timeMults[@]:--1}"; do

    [[ "$betaMult" == "-1" ]] && { bm=""; j_bm=""; } || { bm="inputStruct.betaMult=$betaMult;"; j_bm="-beta=$betaVal"; }
    [[ "$timeMult" == "-1" ]] && { tm=""; j_tm=""; } || { tm="inputStruct.timeMult=$timeMult;"; j_tm="-tmult=$timeMult"; }
    [[ $meshType == "" ]]     && { mt=""; j_mt=""; } || { mt="inputStruct.meshType=\"$meshType\";"; j_mt="-$meshType"; }
    [[ "$meshLvl"  == "-1" ]] && { ml=""; j_ml=""; } || { ml="inputStruct.meshLevel=$meshLvl;"; j_ml="-lev$meshLvl"; }

    baseData=$(sed  -e "s/\[\[PB\]\]/\"$pbName\"/"  \
                    -e "s/\[\[FI\]\]/\"$formID\"/"  \
                    -e "s/\[\[TM\]\]/$tm/"          \
		            -e "s/\[\[BM\]\]/$bm/"          \
		            -e "s/\[\[MT\]\]/$mt/"      \
		            -e "s/\[\[ML\]\]/$ml/"          \
		            <<< "$baseTemplate")

    if [ ${#quadIDs[@]} -ne 0 ]; then
        for quadID in "${quadIDs[@]}"; do
            cmdString="${baseData} inputStruct.quadID=$quadID;"
            cmdBuffer+=("$cmdString")

            jobName="${pbName}-${formID}${j_mt}${j_ml}${j_tm}${j_bm}_${quadID}\'"
            nameBuffer+=("$jobName")
        done
        continue
    fi

    for quadType in "${quadTypes[@]:-""}"; do
    for numSRext in "${numsSRext[@]:--1}"; do
    for numGHext in "${numsGHext[@]:--1}"; do
    for numSRint in "${numsSRint[@]:--1}"; do
    for numGHint in "${numsGHint[@]:--1}"; do
    for numSNGLR in "${numsSNGLR[@]:--1}"; do
    for numBOUND in "${numsBOUND[@]:--1}"; do
    
        [[ $quadType == "" ]] && { qt=""; j_qt=""; } || { qt="inputStruct.quadType=\"$quadType\";"; j_qt="$quadType"; }
        [[ $numSRext -eq -1 ]] && { se=""; j_se=""; } || { se="inputStruct.numSRext=$numSRext;"; j_se="se$numSRext"; }
        [[ $numGHext -eq -1 ]] && { ge=""; j_ge=""; } || { ge="inputStruct.numGHext=$numGHext;"; j_ge="ge$numGHext"; }
        [[ $numSRint -eq -1 ]] && { si=""; j_si=""; } || { si="inputStruct.numSRint=$numSRint;"; j_si="si$numSRint"; }
        [[ $numGHint -eq -1 ]] && { gi=""; j_gi=""; } || { gi="inputStruct.numGHint=$numGHint;"; j_gi="gi$numGHint"; }
        [[ $numSNGLR -eq -1 ]] && { sg=""; j_sg=""; } || { sg="inputStruct.numSNGLR=$numSNGLR;"; j_sg="sg$numSNGLR"; }
        [[ $numBOUND -eq -1 ]] && { bd=""; j_bd=""; } || { bd="inputStruct.numBOUND=$numBOUND;"; j_bd="bd$numBOUND"; }

        quadData=$(sed  -e "s/\[\[QT\]\]/$qt/" \
                        -e "s/\[\[SE\]\]/$se/" \
                        -e "s/\[\[GE\]\]/$ge/" \
                        -e "s/\[\[SI\]\]/$si/" \
                        -e "s/\[\[GI\]\]/$gi/" \
                        -e "s/\[\[SG\]\]/$sg/" \
                        -e "s/\[\[BD\]\]/$bd/" \
                        <<< "$quadTemplate")

        cmdString="${baseData} ${quadData}"
	    cmdBuffer+=("$cmdString")

        jobName="${pbName}-${formID}${j_mt}${j_ml}${j_tm}${j_bm}_${j_qt}${j_se}${j_ge}${j_si}${j_gi}${j_sg}${j_bd}"
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

    if $poolFlag; then
        cmd="${cmd} inputStruct.poolFlag=true;"
    fi
    if $aBldFlag; then
        cmd="${cmd} inputStruct.aBldFlag=true;"
    fi
    if $plotFigs; then
        cmd="${cmd} inputStruct.plotFlag=true;"
    fi
    if $saveFigs; then
        cmd="${cmd} inputStruct.saveFlag=true;"
    fi
    if $saveTemp; then
        cmd="${cmd} inputStruct.saveTemp=true;"
    fi
    echo "Submitting job $name"

    export COMMAND="$cmd"
    sbatch  --job-name="$name"              \
            --output="logs/%j_%x.log"       \
            --qos=gpu                       \
            --partition=gpu                 \
            --nodes=1                       \
            --ntasks-per-node=$poolSize     \
            --mem=$RAMsize                  \
            --gres=gpu:$GPUtype:$GPUcount   \
            --time=1-00:00:00               \
            "./launcher.sh"
done
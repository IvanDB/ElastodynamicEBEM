#!/bin/bash

basePath=$(pwd)
cmdString="basePath=${basePath};"

source "$basePath/parameterLists.config"

# module load matlab

for pbName in "${pbNames[@]:-"def"}"; do

    if [[ $pbName == "def" ]]; then
        echo "Select at least one problem"
        #break
    fi

for meshLvl in "${meshLvls[@]:-1}"; do
for betaVal in "${betaVals[@]:-1}"; do

    pbData="${pbName}_${meshLvl}_${betaVal}"
    cmdString="${cmdString} inputFileName=${pbName}; meshLev=${meshLvl}; beta=${betaVal};"

    for numSRext in "${numsSRext[@]:-16}"; do
    for numGHext in "${numsSRext[@]:-1}"; do
    for numSRint in "${numsSRint[@]:-64}"; do
    for numGHint in "${numsSRint[@]:-3}"; do
    for numSNGLR in "${numsSNGLR[@]:-256}"; do
    for numBOUND in "${numsBOUND[@]:-256}"; do
        
        quadData="${numSRint}-${numGHint}_${numSRext}-${numGHext}_${numSNGLR}_${numBOUND}"
        cmdString="${cmdString} quadString=\"${quadData}\";"

        if $plotFigs; then
            cmdString="${cmdString} plotFlag=true;"
        fi
        if $saveMatx; then
            cmdString="${cmdString} tmpFlag=true;"
        fi
        
        cmdString="${cmdString} main"
        #sbatch --export=COMMAND=${cmdString} launcher.sh

        export COMMAND=${cmdString}
        bash "./launcher.sh"
    done
    done
    done
    done
    done
    done

done
done
done

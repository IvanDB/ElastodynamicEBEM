function problemFileName = constructProblemFileName(pbIndex, inputStruct)
    arguments
        pbIndex (1, 1) {mustBeNonnegative, mustBeInteger} = 0
        inputStruct (1, 1) struct = struct()
    end

    if pbIndex > 0
        problemFileName = getProblemFileName(pbIndex);
        return
    end

    if inputStruct == struct()
        error("Input error", "No input structure provided, and no problem index specified. Please provide at least one of the two inputs.")
    end

    problemFileName = sprintf("input_%s-%s_lev%d.txt", inputStruct.problemName, inputStruct.meshType, inputStruct.meshLevel);
end

function problemFileName = getProblemFileName(pbIndex)
    warning("WIP...")
    switch pbIndex
        case 1
            problemFileName = "input_barH1-symm_lev1.txt";
        case 2
            problemFileName = "input_barH1-symm_lev2.txt";
        case 3
            problemFileName = "input_barH1-symm_lev3.txt";
        otherwise
            error("Input error", "Unrecognized problem index. No corresponding input file name can be generated.")
    end
end
function problemFileName = constructProblemFileName(pbIndex, pbSpecs)
    arguments
        pbIndex  (1, 1) double {mustBeNonnegative, mustBeInteger} = 0
        pbSpecs.pbName (1, 1) string = ""
    end

    pbFileNameList = ["barH1", "barH3", "DesCop-sphere", "DesCop-cube"];

    if pbIndex > 0
        assert(pbIndex <= length(pbFileNameList), sprintf("Unknown problem. ID must be <= %d.", length(pbFileNameList)))
        problemFileName = pbFileNameList(pbIndex) + ".txt";
        return
    end

    if pbSpecs.pbName == ""
        error("No problem specifications provided, and no problem index specified. Please provide at least one of the two inputs.")
    end

    problemFileName = pbSpecs.pbName + ".txt";
end


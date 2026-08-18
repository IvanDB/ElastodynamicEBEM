function problemFileName = constructProblemFileName(pbIndex, pbSpecs)
%CONSTRUCTPROBLEMFILENAME  Resolve the problem (.txt) file name, by preset index or by name.
%   PROBLEMFILENAME = CONSTRUCTPROBLEMFILENAME(PBINDEX, PBSPECS) returns
%   "<name>.txt". If PBINDEX is a positive integer, <name> is looked up in an
%   internal preset list ("barH1", "barH3", "DesCop-sphere", "DesCop-cube").
%   Otherwise, <name> is taken from PBSPECS.pbName, which must then be non-empty.
%
%   Input arguments:
%       PBINDEX - (nonnegative integer, default 0) 1-based index into the
%                 internal preset list; 0 means "use PBSPECS.pbName" instead.
%       PBSPECS - (name-value) pbName (string, default "")
%                 problem name used when PBINDEX is 0.
%
%   Output arguments:
%       PROBLEMFILENAME - (string) file name, without directory,
%                         expected under BASEPATH/inputFiles/.
%
%   Notes:
%       Errors if PBINDEX is 0 and PBSPECS.pbName is empty.
%
%   See also READINPUTFILE, CONSTRUCTMESHFILENAME

arguments (Input)
    pbIndex  (1, 1) double {mustBeNonnegative, mustBeInteger} = 0
    pbSpecs.pbName (1, 1) string = ""
end

arguments (Output)
    problemFileName (1, 1) string
end

pbFileNameList = ["barH1", "barH3", "DesCop-sphere", "DesCop-cube"];

if pbIndex > 0
    assert(pbIndex <= length(pbFileNameList), sprintf("Unknown problem. ID must be <= %d.", length(pbFileNameList)))
    problemFileName = pbFileNameList(pbIndex) + ".txt";
    return
end

assert(pbSpecs.pbName ~= "", "No problem specifications provided, and no problem index specified. Please provide at least one of the two inputs.")
problemFileName = pbSpecs.pbName + ".txt";
end


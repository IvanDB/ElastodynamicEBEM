function meshFileName = constructMeshFileName(pbParam, meshSpecs)
%CONSTRUCTMESHFILENAME  Build the mesh file name for a given problem and refinement level.
%   MESHFILENAME = CONSTRUCTMESHFILENAME(PBPARAM, MESHSPECS) assembles the file name
%   "<meshName>[_<meshType>]_<meshLevel>.mesh" that READSPACEMESH expects, where
%   meshName comes from PBPARAM.meshName. If MESHSPECS.meshType is left empty, the
%   problem's default mesh type (PBPARAM.defaultValues.meshType) is used instead.
%
%   Input arguments:
%       PBPARAM   - (struct) must contain meshName and
%                   defaultValues.meshType, see READINPUTFILE.
%       MESHSPECS - (name-value) meshType (string, default "", meaning "use the
%                   problem default"); meshLevel (nonnegative integer, default 1).
%
%   Output arguments:
%       MESHFILENAME - (string) file name, without directory,
%                      expected under BASEPATH/mesh/<domainName>/.
%
%   See also READSPACEMESH, CONSTRUCTPROBLEMFILENAME

arguments (Input)
    pbParam             (1, 1) struct
    meshSpecs.meshType  (1, 1) string {mustBeText} = "" 
    meshSpecs.meshLevel (1, 1) double {mustBeInteger, mustBeNonnegative} = 1
end

arguments (Output)
    meshFileName (1, 1) string
end

%If no mesh type is provided we use the default one for the problem
if(meshSpecs.meshType == "")
    meshSpecs.meshType = pbParam.defaultValues.meshType;
end

meshFileName = pbParam.meshName;
if(meshSpecs.meshType ~= "")
    meshFileName = meshFileName + "_" + meshSpecs.meshType;
end
meshFileName = meshFileName + "_" + meshSpecs.meshLevel + ".mesh";
end
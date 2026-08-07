function fullFileNames = generateFilenames(basePath, formID, pbParam, domainMesh, quadID)
%GENERATEFILENAMES  Build the on-disk paths used to cache matrices and save the solution.
%   FULLFILENAMES = GENERATEFILENAMES(BASEPATH, FORMID, PBPARAM, DOMAINMESH, QUADID) builds
%   a base file name that encodes the problem (PBPARAM.domainName), the chosen formulation
%   (FORMID), the mesh (DOMAINMESH.name/lev), the time-discretization multipliers
%   (PBPARAM.beta/tMlt) and the quadrature scheme (QUADID), and derives from it the two file
%   paths used by the TIMEMARCHING* routines: one under "tempData" for the intermediate
%   block-Toeplitz matrices/betas, and one under "outputData" for the final solution.
%
%   Input arguments:
%       BASEPATH   - (string) project root.
%       FORMID     - (string) formulation identifier ("ID", "DD", "DN", "DNc" or "IN").
%       PBPARAM    - (struct) physical/time-discretization parameters
%                    (domainName, beta, tMlt), see READINPUTFILE.
%       DOMAINMESH - (struct) triangulated boundary mesh (name, lev), see READSPACEMESH.
%       QUADID     - (string) quadrature scheme identifier, see GENERATEQUADDATA.
%
%   Output arguments:
%       FULLFILENAMES - (struct) with fields tmpFullFilename (under BASEPATH/tempData)
%                       and outFullFilename (under BASEPATH/outputData).
%
%   See also TIMEMARCHINGID, TIMEMARCHINGDD,
%   TIMEMARCHINGDN, TIMEMARCHINGDN_C, TIMEMARCHINGIN
arguments
    basePath    (1, 1) string
    formID      (1, 1) string
    pbParam     (1, 1) struct
    domainMesh  (1, 1) struct
    quadID      (1, 1) string
end


baseFileName = pbParam.domainName + "_" + formID + "_mesh=" + domainMesh.name + domainMesh.lev ...
                        + "_bMult=" + pbParam.beta + "_TMult=" + pbParam.tMlt ...
                        + "_quad=" + quadID;

tmpvarName = "matrix";
outvarName = "solution";
fullFileNames.tmpFullFilename = fullfile(basePath, "tempData", baseFileName + "_" + tmpvarName + ".mat");
fullFileNames.outFullFilename = fullfile(basePath, "outputData", baseFileName + "_" + outvarName + ".mat");

end
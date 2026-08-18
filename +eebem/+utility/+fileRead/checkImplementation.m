function checkImplementation(pbParam)
%CHECKIMPLEMENTATION  Validate that the problem's boundary-integral-equation settings are implemented.
%   CHECKIMPLEMENTATION(PBPARAM) raises an error if PBPARAM.BIE is not "EFIE"
%   (the only formulation family currently implemented), or if PBPARAM.BIE
%   is "EFIE" but PBPARAM.BOU is neither "DIR" nor "NEU".
%   The other BIE values ("MFIE", "HYPE", "ENGI") are recognized as
%   planned-but-not-yet-implemented and reported as such.
%
%   Input arguments:
%       PBPARAM - (struct) must contain the fields BIE and BOU, as read by
%                 READINPUTFILE from the problem file, see READINPUTFILE.
%
%   See also READINPUTFILE

arguments (Input)
    pbParam (1, 1) struct
end

switch pbParam.BIE
    case "EFIE"
        assert(strcmp(pbParam.BOU, "DIR") || strcmp(pbParam.BOU, "NEU"), "Error. Invalid boundary condition for EFIE.")

    case "MFIE"
        assert(false, "Error. MFIE - Method not implemented. Coming Soon.");

    case "HYPE"
        assert(false, "Error. HYPE - Method not implemented. Coming Soon.");

    case "ENGI"
        assert(false, "Error. ENGI - Method not implemented. Coming Soon.");

    otherwise
        assert(false, "Error. Invalid BIE type.");
end

return
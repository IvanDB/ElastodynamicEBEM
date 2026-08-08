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

arguments
    pbParam (1, 1) struct
end

err_flag = 0;
message = "";

switch pbParam.BIE
    case "EFIE"
        if ~strcmp(pbParam.BOU, "DIR") && ~strcmp(pbParam.BOU, "NEU")
            message = "Error. Invalid boundary condition for EFIE.";
            err_flag = 1;
        end
    case "MFIE" 
        message = "Error. MFIE - Method not implemented. Coming Soon.";
        err_flag = 1;
    case "HYPE"
        message = "Error. HYPE - Method not implemented. Coming Soon.";
        err_flag = 1; 
    case "ENGI" 
        message = "Error. ENGI - Method not implemented. Coming Soon.";
        err_flag = 1;     
    otherwise
        message = "Error. Invalid BIE type.";
        err_flag = 1;
end

if err_flag
    error(message)
end
return
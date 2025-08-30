function checkImplementation(pbParam)
% INPUT 
%   - pbParam: struct contenente i parametri del problema

err_flag = 0; 
message = "";

switch pbParam.BIE
    case "EFIE"
        if ~strcmp(pbParam.BOU, "DIR") && ~strcmp(pbParam.BOU, "TEST") && ~strcmp(pbParam.BOU, "NEU")
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
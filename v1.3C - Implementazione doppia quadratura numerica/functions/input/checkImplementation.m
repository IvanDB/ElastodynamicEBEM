function [err_flag, message] = checkImplementation(pb_param)
% INPUT 
%   - pb_param: struct contenente i parametri del problema
%
% OUTPUT:
%   - err_flag: flag intera per gestione errori (default 0)
%   - message: stringa contente eventuale messaggio di errore (default stringa
%   vuota)

err_flag = 0; 
message = '';

switch pb_param.BIE

    case 'EFIE'
        if ~strcmp(pb_param.BOU, 'DIR') && ~strcmp(pb_param.BOU, 'TEST')
            message = 'Error. Invalid boundary condition for EFIE.';
            err_flag = 1;
        end

    case 'MFIE' 
        message = 'Error. MFIE - Method not implemented. Coming Soon.';
        err_flag = 1;

    case 'HYPE'
        message = 'Error. HYPE - Method not implemented. Coming Soon.';
        err_flag = 1; 

    case 'ENGI' 
        message = 'Error. ENGI - Method not implemented. Coming Soon.';
        err_flag = 1;     

    otherwise
        message = 'Error. Invalid BIE type.';
        err_flag = 1;
end

return
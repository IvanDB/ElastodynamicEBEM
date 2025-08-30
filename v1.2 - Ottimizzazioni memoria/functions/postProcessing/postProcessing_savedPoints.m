function [x_val, t_val, i_val] = postProcessing_savedPoints(pb_param)
% INPUT
% - 
% - 
% OUTPUT
% - 

%Selezione del caso di studio per componente x
switch pb_param.domain_type
    case 'screenTest'
        x_val = [];     %Post-processing non presente
        t_val = [];
        i_val = [];
    case 'screenUniform'
        x_val = [];     %Post-processing non presente
        t_val = [];     
        i_val = 3;  
    case 'screenGraded'
        x_val = [];     %Post-processing non presente
        t_val = [];     
        i_val = 3;      
    case 'sphereUniform'
        x_val = [];     %Post-processing non presente
        t_val = [];
        i_val = [];    
    case 'sphereNotUniform'
        x_val = [];     %Post-processing non presente
        t_val = [];
        i_val = [];
    case 'barH1'
        x_val = [0, 0, 0.75;
                 0, 0, 0.50;
                 0, 0, 0.25];
        t_val = linspace(0, pb_param.T_fin, 100*pb_param.T_fin);
        i_val = 3;     
    case 'barH3'
        x_val = [0.25, -0.25, 2];
        t_val = linspace(0, pb_param.T_fin, 100*pb_param.T_fin);
        i_val = 3;     
    case 'sphereWave'
        x_val = [];     %Da implementare se si decide di testare questo caso
        t_val = [];
        i_val = 3;     
    otherwise
        error('Caso ancora non riportato/Errore nei dati')
end


return
end


function pb_param = read_input_file(input_file_name)

% INPUT 
%   - input_file_name: stringa contenente il nome del file di input
%
% OUTPUT:
%   - pb_param: struct contenente i parametri del problema

%% APERTURA FILE
    input_file = fopen(strcat('./input_files/', input_file_name), 'r');
    if ~input_file
        error("Impossibile aprire file di input")
        return
    end

%% PARAMETRI FISICI del PROBLEMA

    %Lettura DENSITA' di MASSA
    line = fgets(input_file);                          %Lettura riga di commento
    pb_param.rho = sscanf(fgets(input_file), '%f');    %Lettura valore
        
    %Lettura MODULO di TAGLIO (mu)
    line = fgets(input_file);                          %Lettura riga di commento
    pb_param.mu = sscanf(fgets(input_file), '%f');     %Lettura valore
    
    %Lettura RAPPORTO di POISSON (nu)
    line = fgets(input_file);                          %Lettura riga di commento
    pb_param.nu = sscanf(fgets(input_file), '%f');     %Lettura valore
        
    %Lettura PARAMETRO di LAME' (lambda)
    line = fgets(input_file);                          %Lettura riga di commento
    pb_param.lambda = sscanf(fgets(input_file), '%f'); %Lettura valore
    %lamba = 2*mu*nu / (1 - 2*nu)
    
    %Calcolo della VELOCITA' c_P delle ONDE P
    pb_param.velP = sqrt((pb_param.lambda+2*pb_param.mu)/pb_param.rho);
    %c_P = sqrt((lamba + 2*mu) / rho)
        
    %Calcolo della VELOCITA' c_S delle ONDE S
    pb_param.velS = sqrt(pb_param.mu/pb_param.rho);
    %c_S = sqrt(mu / rho) 
     
%% PARAMETRI DISCRETIZZAZIONE TEMPORALE
   
    %Lettura ISTANTE TEMPO FINALE                
    line = fgets(input_file);                           %Lettura riga di commento
    pb_param.T_fin = sscanf(fgets(input_file), '%d');   %Lettura valore

    %Lettura NUMERO SOTTOINTERVALLI TEMPORALI
    line = fgets(input_file);                           %Lettura riga di commento
    pb_param.Nt = sscanf(fgets(input_file), '%d');      %Lettura valore
    
%% PARAMETRI DISCRETIZZAZIONE SPAZIALE

    %Lettura del NOME del FILE MESH
    line = fgets(input_file);                               %Lettura riga di commento
    pb_param.domain_type = sscanf(fgets(input_file), '%s'); %Lettura valore
    
    %LIVELLO DI RAFFINAMENTO
    line = fgets(input_file);                           %Lettura riga di commento
    pb_param.lev = sscanf(fgets(input_file), '%d');     %Lettura valore

    %NUMERO NODI GAUSS-HAMMER
    line = fgets(input_file);                           %Lettura riga di commento
    pb_param.ngh = sscanf(fgets(input_file), '%d');     %Lettura valore
%% PARAMETRI EQUAZIONE INTEGRALE

    %Lettura del TIPO di EQUAZIONE INTEGRALE
    line = fgets(input_file);                           %Lettura riga di commento
    pb_param.BIE = sscanf(fgets(input_file), '%s');      %Lettura valore
    
    %Lettura del TIPO di DATO al BORDO
    line = fgets(input_file);                           %Lettura riga di commento
    pb_param.BOU = sscanf(fgets(input_file), '%s');      %Lettura valore
    
    if(pb_param.BOU == "TEST")
	line = fgets(input_file);                           %Lettura riga di commento
    	pb_param.MTDTN = sscanf(fgets(input_file),'%s');    %Lettura valore
    end

 %% CHIUSURA FILE
    fclose(input_file);
return
function pbParam = readInputFile(problemFileName)
% INPUT 
%   - problemFileName: stringa contenente il nome del file di input
%               contenete i dati del problema
%   - methodFileName: stringa contenente il nome del file di input
%               contenete i dati del metodo risolutivo
%
% OUTPUT:
%   - pbParam: struct contenente i parametri del problema

%% APERTURA FILE PROBLEMA
problemFile = fopen(strcat("./inputFiles/", problemFileName), 'r');
if problemFile == -1
    error("Impossibile aprire file coi dati del problema")
end

%% LETTURA PARAMETRI FISICI del PROBLEMA
%Lettura DENSITA' di MASSA
fgets(problemFile);                          %Lettura riga di intestazione
pbParam.rho = sscanf(fgets(problemFile), '%f');    %Lettura valore
    
%Lettura MODULO di TAGLIO (mu)
fgets(problemFile);                          %Lettura riga di intestazione
pbParam.mu = sscanf(fgets(problemFile), '%f');     %Lettura valore

%Lettura RAPPORTO di POISSON (nu)
fgets(problemFile);                          %Lettura riga di intestazione
pbParam.nu = sscanf(fgets(problemFile), '%f');     %Lettura valore
    
%Lettura PARAMETRO di LAME' (lambda)
fgets(problemFile);                          %Lettura riga di intestazione
pbParam.lambda = sscanf(fgets(problemFile), '%f'); %Lettura valore
%lamba = 2*mu*nu / (1 - 2*nu)

%Calcolo della VELOCITA' c_P delle ONDE P
pbParam.velP = sqrt((pbParam.lambda+2*pbParam.mu)/pbParam.rho);
%c_P = sqrt((lamba + 2*mu) / rho)
    
%Calcolo della VELOCITA' c_S delle ONDE S
pbParam.velS = sqrt(pbParam.mu/pbParam.rho);
%c_S = sqrt(mu / rho) 
     
%% LETTURA PARAMETRI DISCRETIZZAZIONE TEMPORALE
%Lettura ISTANTE TEMPO FINALE                
fgets(problemFile);                                 %Lettura riga di intestazione
pbParam.Tfin = sscanf(fgets(problemFile), '%f');   %Lettura valore

%Lettura NUMERO SOTTOINTERVALLI TEMPORALI
fgets(problemFile);                                 %Lettura riga di intestazione
pbParam.nT = sscanf(fgets(problemFile), '%d');      %Lettura valore
    
%% LETTURA PARAMETRI DISCRETIZZAZIONE SPAZIALE
%Lettura del NOME del FILE MESH
fgets(problemFile);                                 %Lettura riga di intestazione
pbParam.domainType = sscanf(fgets(problemFile), '%s'); %Lettura valore

%Lettura LIVELLO RAFFINAMENTO
fgets(problemFile);                                 %Lettura riga di intestazione
pbParam.lev = sscanf(fgets(problemFile), '%d');     %Lettura valore

%% LETTURA PARAMETRI EQUAZIONE INTEGRALE
%Lettura del TIPO di EQUAZIONE INTEGRALE
fgets(problemFile);                                 %Lettura riga di intestazione
pbParam.BIE = sscanf(fgets(problemFile), '%s');     %Lettura valore

%Lettura del TIPO di DATO al BORDO
fgets(problemFile);                           %Lettura riga di intestazione
pbParam.BOU = sscanf(fgets(problemFile), '%s');      %Lettura valore

if(pbParam.BOU == "TEST")
    fgets(problemFile);                           %Lettura riga di intestazione
    pbParam.MTDTN = sscanf(fgets(problemFile),'%s');    %Lettura valore
end

%% CHIUSURA FILE
fclose(problemFile);
return
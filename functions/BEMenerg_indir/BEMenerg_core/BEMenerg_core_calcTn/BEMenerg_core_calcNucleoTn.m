function ris = BEMenerg_core_calcNucleoTn(pbParam, puntoT, t0, t1, ~)
%INPUT
% - ...
% - ...
% OUTPUT
% - ...

%Inizializzazione VETTORE contenente VALUTAZIONE PUNTUALE del DATO al BORDO
ris = zeros(3, 1); 

%Selezione del caso di studio
switch pbParam.BOU
    case 'TEST'
        %Caso test con soluzione phi = (1,1,1) su screen quadrata con L = 1
        L = 1;
        
        %Estrazione VELOCITÀ ONDE S
        velS = pbParam.velS;
        
        %Estrazione VELOCITÀ ONDE P
        velP = pbParam.velP;

        %Calcolo ASCISSA e ORDINATA del NODO di Gauss-Hammer SORGENTE
        x1 = puntoT(1);
        x2 = puntoT(2);

        %Definizione funzioni r*
        r1 = @(y1, y2) y1 - x1;
        r2 = @(y1, y2) y2 - x2;

        %Definizione funzione r 
        r = @(y1, y2) sqrt(r1(y1,y2).^2 + r2(y1,y2).^2); 

        %Definizione delle funzioni J*_t*
        JP_t0 = @(y1, y2) heaviside(t0 - r(y1,y2)/velP) ./ (velP^2 * r(y1,y2));
        JP_t1 = @(y1, y2) heaviside(t1 - r(y1,y2)/velP) ./ (velP^2 * r(y1,y2));
        
        JS_t0 = @(y1, y2) heaviside(t0 - r(y1,y2)/velS) ./ (velS^2 * r(y1,y2));
        JS_t1 = @(y1, y2) heaviside(t1 - r(y1,y2)/velS) ./ (velS^2 * r(y1,y2));
        
        %Definizione delle funzioni I*_t*
        IP_t0 = @(y1, y2) (heaviside(t0 - r(y1,y2)/velP) .* (t0 - r(y1,y2)/velP) .* (t0 + r(y1,y2)/velP)) ./ (2 * r(y1,y2).^3);
        IP_t1 = @(y1, y2) (heaviside(t1 - r(y1,y2)/velP) .* (t1 - r(y1,y2)/velP) .* (t1 + r(y1,y2)/velP)) ./ (2 * r(y1,y2).^3);
        
        IS_t0 = @(y1, y2) (heaviside(t0 - r(y1,y2)/velS) .* (t0 - r(y1,y2)/velS) .* (t0 + r(y1,y2)/velS)) ./ (2 * r(y1,y2).^3);
        IS_t1 = @(y1, y2) (heaviside(t1 - r(y1,y2)/velS) .* (t1 - r(y1,y2)/velS) .* (t1 + r(y1,y2)/velS)) ./ (2 * r(y1,y2).^3);
        
        %Definizione delle funzioni J*_t*_** 
        JP_t0_11 = @(y1,y2) (r1(y1,y2).^2) .* JP_t0(y1,y2) ./ (r(y1,y2).^2);
        JP_t1_11 = @(y1,y2) (r1(y1,y2).^2) .* JP_t1(y1,y2) ./ (r(y1,y2).^2);
        
        JS_t0_11 = @(y1,y2) (r1(y1,y2).^2) .* JS_t0(y1,y2) ./ (r(y1,y2).^2);
        JS_t1_11 = @(y1,y2) (r1(y1,y2).^2) .* JS_t1(y1,y2) ./ (r(y1,y2).^2);
        
        JP_t0_22 = @(y1,y2) (r2(y1,y2).^2) .* JP_t0(y1,y2) ./ (r(y1,y2).^2);
        JP_t1_22 = @(y1,y2) (r2(y1,y2).^2) .* JP_t1(y1,y2) ./ (r(y1,y2).^2);
        
        JS_t0_22 = @(y1,y2) (r2(y1,y2).^2) .* JS_t0(y1,y2) ./ (r(y1,y2).^2);
        JS_t1_22 = @(y1,y2) (r2(y1,y2).^2) .* JS_t1(y1,y2) ./ (r(y1,y2).^2);
     
        JP_t0_12 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)) .* JP_t0(y1,y2) ./ (r(y1,y2).^2);
        JP_t1_12 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)) .* JP_t1(y1,y2) ./ (r(y1,y2).^2);
        
        JS_t0_12 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)) .* JS_t0(y1,y2) ./ (r(y1,y2).^2);
        JS_t1_12 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)) .* JS_t1(y1,y2) ./ (r(y1,y2).^2);   
        %----------------------------------------------------------------

        %Definizione delle funzioni I*_t*_**   
        IP_t0_11 = @(y1,y2) (r1(y1,y2).^2) .* IP_t0(y1,y2) ./ (r(y1,y2).^2);
        IP_t1_11 = @(y1,y2) (r1(y1,y2).^2) .* IP_t1(y1,y2) ./ (r(y1,y2).^2);
        
        IS_t0_11 = @(y1,y2) (r1(y1,y2).^2) .* IS_t0(y1,y2) ./ (r(y1,y2).^2);
        IS_t1_11 = @(y1,y2) (r1(y1,y2).^2) .* IS_t1(y1,y2) ./ (r(y1,y2).^2);
        
        IP_t0_22 = @(y1,y2) (r2(y1,y2).^2) .* IP_t0(y1,y2) ./ (r(y1,y2).^2);
        IP_t1_22 = @(y1,y2) (r2(y1,y2).^2) .* IP_t1(y1,y2) ./ (r(y1,y2).^2);
        
        IS_t0_22 = @(y1,y2) (r2(y1,y2).^2) .* IS_t0(y1,y2) ./ (r(y1,y2).^2);
        IS_t1_22 = @(y1,y2) (r2(y1,y2).^2) .* IS_t1(y1,y2) ./ (r(y1,y2).^2);
        
        IP_t0_12 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)) .* IP_t0(y1,y2) ./ (r(y1,y2).^2);
        IP_t1_12 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)) .* IP_t1(y1,y2) ./ (r(y1,y2).^2);
        
        IS_t0_12 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)) .* IS_t0(y1,y2) ./ (r(y1,y2).^2);
        IS_t1_12 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)) .* IS_t1(y1,y2) ./ (r(y1,y2).^2);
        %----------------------------------------------------------------
         
         
        %Definizione dei nuclei integrali  
        
        nucl_1 = @(y1,y2) - (JP_t0_11(y1,y2) - JS_t0_11(y1,y2)) ...
                            + (IP_t0(y1,y2) - IS_t0(y1,y2)) ...
                            - 3*(IP_t0_11(y1,y2) - IS_t0_11(y1,y2)) ...
                            - JS_t0(y1,y2) ...
                          + (JP_t1_11(y1,y2) - JS_t1_11(y1,y2)) ...
                            - (IP_t1(y1,y2) - IS_t1(y1,y2)) ...
                            + 3*(IP_t1_11(y1,y2) - IS_t1_11(y1,y2)) ...
                            + JS_t1(y1,y2) ...
                          - (JP_t0_12(y1,y2) - JS_t0_12(y1,y2)) ...
                            - 3*(IP_t0_12(y1,y2) - IS_t0_12(y1,y2)) ...
                          + (JP_t1_12(y1,y2) - JS_t1_12(y1,y2)) ...
                            + 3*(IP_t1_12(y1,y2) - IS_t1_12(y1,y2));

        
        nucl_2 = @(y1,y2) - (JP_t0_12(y1,y2) - JS_t0_12(y1,y2)) ...
                            - 3*(IP_t0_12(y1,y2) - IS_t0_12(y1,y2)) ...
                          + (JP_t1_12(y1,y2) - JS_t1_12(y1,y2)) ...
                            + 3*(IP_t1_12(y1,y2) - IS_t1_12(y1,y2)) ...
                          - (JP_t0_22(y1,y2) - JS_t0_22(y1,y2)) ...
                            + (IP_t0(y1,y2) - IS_t0(y1,y2)) ...
                            - 3*(IP_t0_22(y1,y2) - IS_t0_22(y1,y2)) ...
                            - JS_t0(y1,y2)...
                          + (JP_t1_22(y1,y2) - JS_t1_22(y1,y2)) ...
                            - (IP_t1(y1,y2) - IS_t1(y1,y2)) ...
                            + 3*(IP_t1_22(y1,y2) - IS_t1_22(y1,y2)) ...
                            + JS_t1(y1,y2);
                     
        nucl_3 = @(y1,y2) + (IP_t0(y1,y2) - IS_t0(y1,y2)) - JS_t0(y1,y2) ...
                          - (IP_t1(y1,y2) - IS_t1(y1,y2)) + JS_t1(y1,y2);
        
        %Calcolo integrali
        switch pbParam.MTDTN
            case 'classic'
    	        ris = BEMenerg_core_TnTest_classic(nucl_1, nucl_2, nucl_3, L);
            case 'doppioIntergral1D'
    	        ris = BEMenerg_core_TnTest_doppioIntegral1D(nucl_1, nucl_2, nucl_3, L);
            case 'doppioIntergral1DSpezzato'
	            ris = BEMenerg_core_TnTest_doppioIntegral1DSpezzato(nucl_1, nucl_2, nucl_3, L, x1, x2);
        end
	
        %Applicazione coeffciente comune
        ris = ris ./ (4*pi*pbParam.rho);
    case 'DIR'
        %Selezione dato al bordo
        switch pbParam.domainType
            case 'sphereUniform'
                g = @(x, t) [0; 0; 0.5 .* t]; 

            case 'sphereNotUniform'
                g = @(x, t) [0; 0; 0.5 .* t]; 

            case 'screenUniform'
                g = @(x, t) [0; 0; ([sin(t).^5, 1] * [t < pi/2; t >= pi/2]) .* (x(1) .* x(1))];

            case 'screenGraded'
                g = @(x, t) [0; 0; ([sin(t).^5, 1] * [t < pi/2; t >= pi/2]) .* (x(1) .* x(1))];

            case 'barH1'
                h = 1;
                velP = pbParam.velP;
                hat_k = ceil((velP * pbParam.Tfin) / (2 .* h)) - 1;
                k_val = (0 : hat_k)';
                sol_an = @(x, t) sum((-1).^k_val .* ( (velP.*t - 2.*h.*k_val - (h-x(3))) .* ((velP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                                - (velP.*t - 2.*h.*(k_val+1) + (h-x(3))) .* ((velP.*t - 2.*h.*(k_val+1) + (h-x(3))) > 0)));
                g = @(x, t) [0; 0; sol_an(x, t)];

            case 'barH3'
                h = 3;
                velP = pbParam.velP;
                hat_k = ceil((velP * pbParam.Tfin) / (2 .* h)) - 1;
                k_val = (0 : hat_k)';
                sol_an = @(x, t) sum((-1).^k_val .* ( (velP.*t - 2.*h.*k_val - (h-x(3))) .* ((velP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                                - (velP.*t - 2.*h.*(k_val+1) + (h-x(3))) .* ((velP.*t - 2.*h.*(k_val+1) + (h-x(3))) > 0)));
                g = @(x, t) [0; 0; sol_an(x, t)]; 

            case 'sphereWave'
                velP = pbParam.velP;
                ondaIncid = @(x, t) exp(-20 .* ((x(1) - 2 + velP.*t - 0.475).^2));
                g = @(x, t) [-ondaIncid(x, t); 0; 0];

            case 'elementoIndustriale'
                velP = pbParam.velP;
                ondaIncid = @(x, t) exp(-20 .* ((x(3) + 1 - velP.*t).^2));
                g = @(x, t) [0; 0; -ondaIncid(x, t)];
                
            otherwise
                error('Caso ancora non riportato/Errore nei dati')
        end
        
        %Calcolo dato al bordo corrente
        ris = g(puntoT, t1) - g(puntoT, t0);
        
    otherwise
        error('Caso ancora non riportato/Errore nei dati')
end
return
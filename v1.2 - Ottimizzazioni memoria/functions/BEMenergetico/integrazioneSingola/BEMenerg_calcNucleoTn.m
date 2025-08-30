function ris = BEMenerg_calcNucleoTn(pb_param, ghS, t0, t1, ind_RHS)
%INPUT
% - ...
% - ...
% OUTPUT
% - ...

%Inizializzazione VETTORE contenente VALUTAZIONE PUNTUALE del DATO al BORDO
ris = zeros(3, 1); 

%Selezione del caso di studio
switch pb_param.BOU
    case 'TEST'
        %Caso test con soluzione phi = (1,1,1) su screen quadrata con L = 1
        L = 1;
        
        %Estrazione VELOCITÀ ONDE S
        velC_S = pb_param.velS;
        
        %Estrazione VELOCITÀ ONDE P
        velC_P = pb_param.velP;

        %Calcolo ASCISSA e ORDINATA del NODO di Gauss-Hammer SORGENTE
        x1 = ghS(1);
        x2 = ghS(2);

        %Definizione funzioni r*
        r1 = @(y1, y2) y1 - x1;
        r2 = @(y1, y2) y2 - x2;

        %Definizione funzione r 
        r = @(y1, y2) sqrt(r1(y1,y2).^2 + r2(y1,y2).^2); 

        %Definizione delle funzioni J*_t*
        JP_t0 = @(y1, y2) heaviside(t0 - r(y1,y2)/velC_P) ./ (velC_P^2 * r(y1,y2));
        JP_t1 = @(y1, y2) heaviside(t1 - r(y1,y2)/velC_P) ./ (velC_P^2 * r(y1,y2));
        
        JS_t0 = @(y1, y2) heaviside(t0 - r(y1,y2)/velC_S) ./ (velC_S^2 * r(y1,y2));
        JS_t1 = @(y1, y2) heaviside(t1 - r(y1,y2)/velC_S) ./ (velC_S^2 * r(y1,y2));
        
        %Definizione delle funzioni I*_t*
        IP_t0 = @(y1, y2) (heaviside(t0 - r(y1,y2)/velC_P) .* (t0 - r(y1,y2)/velC_P) .* (t0 + r(y1,y2)/velC_P)) ./ (2 * r(y1,y2).^3);
        IP_t1 = @(y1, y2) (heaviside(t1 - r(y1,y2)/velC_P) .* (t1 - r(y1,y2)/velC_P) .* (t1 + r(y1,y2)/velC_P)) ./ (2 * r(y1,y2).^3);
        
        IS_t0 = @(y1, y2) (heaviside(t0 - r(y1,y2)/velC_S) .* (t0 - r(y1,y2)/velC_S) .* (t0 + r(y1,y2)/velC_S)) ./ (2 * r(y1,y2).^3);
        IS_t1 = @(y1, y2) (heaviside(t1 - r(y1,y2)/velC_S) .* (t1 - r(y1,y2)/velC_S) .* (t1 + r(y1,y2)/velC_S)) ./ (2 * r(y1,y2).^3);
        
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
	switch pb_param.MTDTN
	    case 'classic'
        	ris = BEMenerg_TnTest_classic(nucl_1, nucl_2, nucl_3, L);
	    case 'doppioIntergral1D'
        	ris = BEMenerg_TnTest_doppioIntegral1D(nucl_1, nucl_2, nucl_3, L);
	    case 'doppioIntergral1DSpezzato'
		    ris = BEMenerg_TnTest_doppioIntegral1DSpezzato(nucl_1, nucl_2, nucl_3, L, x1, x2);
    end
	
        %Applicazione coeffciente comune
        ris = ris ./ (4*pi*pb_param.rho);
    case 'DIR'
        %Selezione dato al bordo
        switch pb_param.domain_type
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
                cP = pb_param.velP;
                hat_k = ceil((pb_param.velP * pb_param.T_fin) / (2 .* h)) - 1;
                k_val = (0 : hat_k)';
                sol_an = @(x, t) sum((-1).^k_val .* ( (cP.*t - 2.*h.*k_val - (h-x(3))) .* ((cP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                                - (cP.*t - 2.*h.*(k_val+1) + (h-x(3))) .* ((cP.*t - 2.*h.*(k_val+1) + (h-x(3))) > 0)));
                
                %Coefficiente moltiplicativo 4 = cP^2 necessario
                % per ottenere coerenza con risultati tesi
                g = @(x, t) [0; 0; sol_an(x, t)];              %dato teorico (e corretto)
                %g = @(x, t) (cP.^2) .* [0; 0; sol_an(x, t)];    %dato modificato
            case 'barH3'
                h = 3;
                cP = pb_param.velP;
                hat_k = ceil((pb_param.velP * pb_param.T_fin) / (2 .* h)) - 1;
                k_val = (0 : hat_k)';
                sol_an = @(x, t) sum((-1).^k_val .* ( (cP.*t - 2.*h.*k_val - (h-x(3))) .* ((cP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                                - (cP.*t - 2.*h.*(k_val+1) + (h-x(3))) .* ((cP.*t - 2.*h.*(k_val+1) + (h-x(3))) > 0)));
                g = @(x, t) [0; 0; sol_an(x, t)]; 
            case 'sphereWave'
                cP = pb_param.velP;
                onda_incid = @(x, t) exp(-20 .* ((x(1) - 2 + cP.*t - 0.475).^2));
                g = @(x, t) [-onda_incid(x, t); 0; 0];
            otherwise
                error('Caso ancora non riportato/Errore nei dati')
        end
        
        %Calcolo dato al bordo corrente
        ris = g(ghS, t1) - g(ghS, t0);
    otherwise
        error('Caso ancora non riportato/Errore nei dati')
end
return
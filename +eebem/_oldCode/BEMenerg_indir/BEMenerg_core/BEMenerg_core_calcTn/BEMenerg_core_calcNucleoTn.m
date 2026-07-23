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
               
            case "DesCop-cube"
                a = 0.1;
                b = 100;
                cP = pbParam.velP;
                cS = pbParam.velS;
                rho = pbParam.rho;
        
                qP = @(r, t) sqrt(a) * (a*b - t + (r / cP));
                qS = @(r, t) sqrt(a) * (a*b - t + (r / cS));
                F = @(r, t) (1 / (2 * a * (r^2))) * (exp(-(qP(r, t)^2)) - exp(-(qS(r, t)^2)) + sqrt(a * pi) * (a*b - t) * (erf(qP(r, t)) - erf(qS(r, t))));
                
                fP = @(r, t) exp(-a * ((t - (r / cP) - (a*b))^2));
                fS = @(r, t) exp(-a * ((t - (r / cS) - (a*b))^2));
        
                d = @(i, j) (i == j);
                x0 = [1, 1, 1];
        
                u_ij = @(r, t, i, j) ((3 * r(i) * r(j) / (norm(r)^3)) - (d(i, j) / norm(r))) * F(norm(r), t) ...
                                    + (r(i) * r(j) / (norm(r)^3)) * ((fP(norm(r), t) / (cP^2)) - (fS(norm(r), t) / (cS^2))) ...
                                    + (d(i, j) / norm(r)) * (fS(norm(r), t) / (cS^2));
        
                u_i = @(r, t, i) u_ij(r, t, i, 1) + u_ij(r, t, i, 2) + u_ij(r, t, i, 3);
                
                u = @(r, t) (1 / (4*pi*rho)) .* [u_i(r, t, 1); u_i(r, t, 2); u_i(r, t, 3)];
        
                g = @(x, t) u(x - x0, t);
                
            case "DesCop-sphere"
                nu = pbParam.nu;
                cP = pbParam.velP;
                rho = pbParam.rho;
                R = 1;
                p0 = 143.5;
                a = 1279;
                b = 12792;
        
                omega0 = (cP * sqrt(1 - 2*nu)) / (R * (1 - nu));
                alpha0 = (cP * (1 - 2*nu)) / (R * (1 - nu));
                A2 = @(z) omega0 + (alpha0 - z)^2;
                A = @(z) sqrt(A2(z));
                beta0 = @(z) asin((alpha0 - z) / A(z));
        
               
                g_func = @(t, z) (p0 / (rho * A2(z) * omega0)) * ...
                    ((1 / cP) * (-z * omega0 * exp(-z*t) + alpha0 * A(z) * exp(-alpha0*t) * cos(omega0*t - beta0(z)) + omega0 *  A(z) * exp(-alpha0*t) * sin(omega0*t - beta0(z))) ...
                      - (1 / R) * (A(z) * exp(-alpha0*t) * cos(omega0*t - beta0(z)) - omega0 * exp(-z*t)));
        
                g = @(x, t) (g_func(t, a) - g_func(t, b)) .* (-x' ./ norm(x)); 
            otherwise
                error('Caso ancora non riportato/Errore nei dati')
        end
        
        %Calcolo dato al bordo corrente
        ris = g(puntoT, t1) - g(puntoT, t0);
        
    otherwise
        error('Caso ancora non riportato/Errore nei dati')
end
return
function ris = time3D_intLE_rhs_screen_y_0(L,pb_param,x1,x3,t0,t1,ind_RHS)

%La function time3D_intLE() prende in INPUT:
%- la variabile L contenente la lunghezza del lato della screen quadrata
%
%- la struct pb_param contenente i parametri del problema
%
%- la variabile x1 contenente l'ascissa del punto sorgente sp
%
%- la variabile x2 contenente l'ordinata del punto sorgente sp
%
%- la variabile t0 che contiene il primo istante temporale
%  (istante di tempo precedente all'istante corrente del time-marching)
%
%- la variabile t1 che contiene il secondo istante temporale 
%  (istante di tempo corrente del time-marching)
%
%- la variabile ind_RHS l'indice che fornisce il tipo di dato al bordo
%  assegnato al triangolo sorgente corrente
%
%e restituisce in OUTPUT:
%- il vettore 3x1 ris contenente il risultato della valutazione puntuale 
%  del dato al bordo
%
%-------------------------------------------------------------------------

%Inizializzazione del VETTORE contenente la VALUTAZIONE PUNTUALE del 
%DATO al BORDO
ris = zeros(3,1); 

%VELOCITÀ delle ONDE S
velC_S = pb_param.velS;

%VELOCITÀ delle ONDE P
velC_P = pb_param.velP;

%Definiamo le COMPONENTI del VETTORE DISTANZA r tra un PUNTO di CAMPO
%e un PUNTO SORGENTE fissato (nodo di Gauss-Hammer)
r1 = @(y1,y3) y1-x1;
r3 = @(y1,y3) y3-x3;

%Definiamo il MODULO del VETTORE DISTANZA r tra un PUNTO di CAMPO
%e il PUNTO SORGENTE fissato (nodo di Gauss-Hammer)
r = @(y1,y3) sqrt(r1(y1,y3).^2+r3(y1,y3).^2);

switch ind_RHS
    
    case 1 
        %Dato che permette di ottenere come soluzione phi = (1,0,0) sulla
        %screen quadrata
        
        %----------------------------------------------------------------
        %Definizione delle funzioni I_*
        funI_Pt0 = @(y1,y3) (heaviside(t0-r(y1,y3)/velC_P).*(t0-r(y1,y3)/velC_P).*(t0+r(y1,y3)/velC_P))./(2*r(y1,y3).^3);
        funI_Pt1 = @(y1,y3) (heaviside(t1-r(y1,y3)/velC_P).*(t1-r(y1,y3)/velC_P).*(t1+r(y1,y3)/velC_P))./(2*r(y1,y3).^3);
        
        funI_St0 = @(y1,y3) (heaviside(t0-r(y1,y3)/velC_S).*(t0-r(y1,y3)/velC_S).*(t0+r(y1,y3)/velC_S))./(2*r(y1,y3).^3);
        funI_St1 = @(y1,y3) (heaviside(t1-r(y1,y3)/velC_S).*(t1-r(y1,y3)/velC_S).*(t1+r(y1,y3)/velC_S))./(2*r(y1,y3).^3);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni J_*
        funJ_St0 = @(y1,y3) heaviside(t0-r(y1,y3)/velC_S)./(velC_S^2*r(y1,y3));
        funJ_St1 = @(y1,y3) heaviside(t1-r(y1,y3)/velC_S)./(velC_S^2*r(y1,y3));
        
        funJ_Pt0 = @(y1,y3) heaviside(t0-r(y1,y3)/velC_P)./(velC_P^2*r(y1,y3));
        funJ_Pt1 = @(y1,y3) heaviside(t1-r(y1,y3)/velC_P)./(velC_P^2*r(y1,y3));
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_I_*
        
        funtI_Pt0_1_1 = @(y1,y3) (r1(y1,y3).^2).*funI_Pt0(y1,y3)./(r(y1,y3).^2);
        funtI_Pt1_1_1 = @(y1,y3) (r1(y1,y3).^2).*funI_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtI_St0_1_1 = @(y1,y3) (r1(y1,y3).^2).*funI_St0(y1,y3)./(r(y1,y3).^2);
        funtI_St1_1_1 = @(y1,y3) (r1(y1,y3).^2).*funI_St1(y1,y3)./(r(y1,y3).^2);
        
        funtI_Pt0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_Pt0(y1,y3)./(r(y1,y3).^2);
        funtI_Pt1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtI_St0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_St0(y1,y3)./(r(y1,y3).^2);
        funtI_St1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_St1(y1,y3)./(r(y1,y3).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_J_*
        
        funtJ_Pt0_1_1 = @(y1,y3) (r1(y1,y3).^2).*funJ_Pt0(y1,y3)./(r(y1,y3).^2);
        funtJ_Pt1_1_1 = @(y1,y3) (r1(y1,y3).^2).*funJ_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtJ_St0_1_1 = @(y1,y3) (r1(y1,y3).^2).*funJ_St0(y1,y3)./(r(y1,y3).^2);
        funtJ_St1_1_1 = @(y1,y3) (r1(y1,y3).^2).*funJ_St1(y1,y3)./(r(y1,y3).^2);
        
        funtJ_Pt0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_Pt0(y1,y3)./(r(y1,y3).^2);
        funtJ_Pt1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtJ_St0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_St0(y1,y3)./(r(y1,y3).^2);
        funtJ_St1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_St1(y1,y3)./(r(y1,y3).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni integrande
        g_1 = @(y1,y3) -(funtJ_Pt0_1_1(y1,y3)-funtJ_St0_1_1(y1,y3))+(funI_Pt0(y1,y3)-funI_St0(y1,y3))-3*(funtI_Pt0_1_1(y1,y3)-funtI_St0_1_1(y1,y3))-funJ_St0(y1,y3)...
            +(funtJ_Pt1_1_1(y1,y3)-funtJ_St1_1_1(y1,y3))-(funI_Pt1(y1,y3)-funI_St1(y1,y3))+3*(funtI_Pt1_1_1(y1,y3)-funtI_St1_1_1(y1,y3))+funJ_St1(y1,y3);

        g_3 = @(y1,y3) -(funtJ_Pt0_1_3(y1,y3)-funtJ_St0_1_3(y1,y3))-3*(funtI_Pt0_1_3(y1,y3)-funtI_St0_1_3(y1,y3))...
            +(funtJ_Pt1_1_3(y1,y3)-funtJ_St1_1_3(y1,y3))+3*(funtI_Pt1_1_3(y1,y3)-funtI_St1_1_3(y1,y3));
        %----------------------------------------------------------------
      
        %Calcolo degli integrali
        ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        ris(3) = integral2(g_3,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
        %----------------------------------------------------------------
        
        ris(1) = ris(1)/(4*pi*pb_param.rho);
        ris(3) = ris(3)/(4*pi*pb_param.rho);
        
    case 2
        %Dato che permette di ottenere come soluzione phi = (0,1,0) sulla
        %screen quadrata
        
        %----------------------------------------------------------------
        %Definizione delle funzioni I_*
        funI_Pt0 = @(y1,y3) (heaviside(t0-r(y1,y3)/velC_P).*(t0-r(y1,y3)/velC_P).*(t0+r(y1,y3)/velC_P))./(2*r(y1,y3).^3);
        funI_Pt1 = @(y1,y3) (heaviside(t1-r(y1,y3)/velC_P).*(t1-r(y1,y3)/velC_P).*(t1+r(y1,y3)/velC_P))./(2*r(y1,y3).^3);
        
        funI_St0 = @(y1,y3) (heaviside(t0-r(y1,y3)/velC_S).*(t0-r(y1,y3)/velC_S).*(t0+r(y1,y3)/velC_S))./(2*r(y1,y3).^3);
        funI_St1 = @(y1,y3) (heaviside(t1-r(y1,y3)/velC_S).*(t1-r(y1,y3)/velC_S).*(t1+r(y1,y3)/velC_S))./(2*r(y1,y3).^3);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni J_*
        funJ_St0 = @(y1,y3) heaviside(t0-r(y1,y3)/velC_S)./(velC_S^2*r(y1,y3));
        funJ_St1 = @(y1,y3) heaviside(t1-r(y1,y3)/velC_S)./(velC_S^2*r(y1,y3));
        
        %----------------------------------------------------------------
        %Definizione delle funzioni integrande
        g_2 = @(y1,y3) +(funI_Pt0(y1,y3)-funI_St0(y1,y3))-funJ_St0(y1,y3)...
            -(funI_Pt1(y1,y3)-funI_St1(y1,y3))+funJ_St1(y1,y3);
        %----------------------------------------------------------------
        %Calcolo degli integrali    
        ris(2) = integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        %----------------------------------------------------------------
        
        
        ris(2) = ris(2)/(4*pi*pb_param.rho);
         
    case 3 
        %Dato che permette di ottenere come soluzione phi = (0,0,1) sulla
        %screen quadrata
        
        %----------------------------------------------------------------
        %Definizione delle funzioni I_*
        funI_Pt0 = @(y1,y3) (heaviside(t0-r(y1,y3)/velC_P).*(t0-r(y1,y3)/velC_P).*(t0+r(y1,y3)/velC_P))./(2*r(y1,y3).^3);
        funI_Pt1 = @(y1,y3) (heaviside(t1-r(y1,y3)/velC_P).*(t1-r(y1,y3)/velC_P).*(t1+r(y1,y3)/velC_P))./(2*r(y1,y3).^3);
        
        funI_St0 = @(y1,y3) (heaviside(t0-r(y1,y3)/velC_S).*(t0-r(y1,y3)/velC_S).*(t0+r(y1,y3)/velC_S))./(2*r(y1,y3).^3);
        funI_St1 = @(y1,y3) (heaviside(t1-r(y1,y3)/velC_S).*(t1-r(y1,y3)/velC_S).*(t1+r(y1,y3)/velC_S))./(2*r(y1,y3).^3);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni J_*
        funJ_St0 = @(y1,y3) heaviside(t0-r(y1,y3)/velC_S)./(velC_S^2*r(y1,y3));
        funJ_St1 = @(y1,y3) heaviside(t1-r(y1,y3)/velC_S)./(velC_S^2*r(y1,y3));
        
        funJ_Pt0 = @(y1,y3) heaviside(t0-r(y1,y3)/velC_P)./(velC_P^2*r(y1,y3));
        funJ_Pt1 = @(y1,y3) heaviside(t1-r(y1,y3)/velC_P)./(velC_P^2*r(y1,y3));
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_I_*
        
        funtI_Pt0_3_3 = @(y1,y3) (r3(y1,y3).^2).*funI_Pt0(y1,y3)./(r(y1,y3).^2);
        funtI_Pt1_3_3 = @(y1,y3) (r3(y1,y3).^2).*funI_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtI_St0_3_3 = @(y1,y3) (r3(y1,y3).^2).*funI_St0(y1,y3)./(r(y1,y3).^2);
        funtI_St1_3_3 = @(y1,y3) (r3(y1,y3).^2).*funI_St1(y1,y3)./(r(y1,y3).^2);
        
        funtI_Pt0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_Pt0(y1,y3)./(r(y1,y3).^2);
        funtI_Pt1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtI_St0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_St0(y1,y3)./(r(y1,y3).^2);
        funtI_St1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_St1(y1,y3)./(r(y1,y3).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_J_*
        
        funtJ_Pt0_3_3 = @(y1,y3) (r3(y1,y3).^2).*funJ_Pt0(y1,y3)./(r(y1,y3).^2);
        funtJ_Pt1_3_3 = @(y1,y3) (r3(y1,y3).^2).*funJ_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtJ_St0_3_3 = @(y1,y3) (r3(y1,y3).^2).*funJ_St0(y1,y3)./(r(y1,y3).^2);
        funtJ_St1_3_3 = @(y1,y3) (r3(y1,y3).^2).*funJ_St1(y1,y3)./(r(y1,y3).^2);
        
        funtJ_Pt0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_Pt0(y1,y3)./(r(y1,y3).^2);
        funtJ_Pt1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtJ_St0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_St0(y1,y3)./(r(y1,y3).^2);
        funtJ_St1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_St1(y1,y3)./(r(y1,y3).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni integrande
        g_1 = @(y1,y3) -(funtJ_Pt0_1_3(y1,y3)-funtJ_St0_1_3(y1,y3))-3*(funtI_Pt0_1_3(y1,y3)-funtI_St0_1_3(y1,y3))...
            +(funtJ_Pt1_1_3(y1,y3)-funtJ_St1_1_3(y1,y3))+3*(funtI_Pt1_1_3(y1,y3)-funtI_St1_1_3(y1,y3));
        
        g_3 = @(y1,y3) -(funtJ_Pt0_3_3(y1,y3)-funtJ_St0_3_3(y1,y3))+(funI_Pt0(y1,y3)-funI_St0(y1,y3))-3*(funtI_Pt0_3_3(y1,y3)-funtI_St0_3_3(y1,y3))-funJ_St0(y1,y3)...
            +(funtJ_Pt1_3_3(y1,y3)-funtJ_St1_3_3(y1,y3))-(funI_Pt1(y1,y3)-funI_St1(y1,y3))+3*(funtI_Pt1_3_3(y1,y3)-funtI_St1_3_3(y1,y3))+funJ_St1(y1,y3);
        %----------------------------------------------------------------
      
        %Calcolo degli integrali
        ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
        ris(3) = integral2(g_3,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        %----------------------------------------------------------------
        
        ris(1) = ris(1)/(4*pi*pb_param.rho);
        ris(3) = ris(3)/(4*pi*pb_param.rho);
        
    case 4
        %Dato che permette di ottenere come soluzione phi = (1,1,1) sulla
        %screen quadrata
        
        %----------------------------------------------------------------
        %Definizione delle funzioni I_*
        funI_Pt0 = @(y1,y3) (heaviside(t0-r(y1,y3)/velC_P).*(t0-r(y1,y3)/velC_P).*(t0+r(y1,y3)/velC_P))./(2*r(y1,y3).^3);
        funI_Pt1 = @(y1,y3) (heaviside(t1-r(y1,y3)/velC_P).*(t1-r(y1,y3)/velC_P).*(t1+r(y1,y3)/velC_P))./(2*r(y1,y3).^3);
        
        funI_St0 = @(y1,y3) (heaviside(t0-r(y1,y3)/velC_S).*(t0-r(y1,y3)/velC_S).*(t0+r(y1,y3)/velC_S))./(2*r(y1,y3).^3);
        funI_St1 = @(y1,y3) (heaviside(t1-r(y1,y3)/velC_S).*(t1-r(y1,y3)/velC_S).*(t1+r(y1,y3)/velC_S))./(2*r(y1,y3).^3);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni J_*
        funJ_St0 = @(y1,y3) heaviside(t0-r(y1,y3)/velC_S)./(velC_S^2*r(y1,y3));
        funJ_St1 = @(y1,y3) heaviside(t1-r(y1,y3)/velC_S)./(velC_S^2*r(y1,y3));
        
        funJ_Pt0 = @(y1,y3) heaviside(t0-r(y1,y3)/velC_P)./(velC_P^2*r(y1,y3));
        funJ_Pt1 = @(y1,y3) heaviside(t1-r(y1,y3)/velC_P)./(velC_P^2*r(y1,y3));
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_I_*   
        funtI_Pt0_1_1 = @(y1,y3) (r1(y1,y3).^2).*funI_Pt0(y1,y3)./(r(y1,y3).^2);
        funtI_Pt1_1_1 = @(y1,y3) (r1(y1,y3).^2).*funI_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtI_St0_1_1 = @(y1,y3) (r1(y1,y3).^2).*funI_St0(y1,y3)./(r(y1,y3).^2);
        funtI_St1_1_1 = @(y1,y3) (r1(y1,y3).^2).*funI_St1(y1,y3)./(r(y1,y3).^2);
        
        funtI_Pt0_3_3 = @(y1,y3) (r3(y1,y3).^2).*funI_Pt0(y1,y3)./(r(y1,y3).^2);
        funtI_Pt1_3_3 = @(y1,y3) (r3(y1,y3).^2).*funI_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtI_St0_3_3 = @(y1,y3) (r3(y1,y3).^2).*funI_St0(y1,y3)./(r(y1,y3).^2);
        funtI_St1_3_3 = @(y1,y3) (r3(y1,y3).^2).*funI_St1(y1,y3)./(r(y1,y3).^2);
        
        funtI_Pt0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_Pt0(y1,y3)./(r(y1,y3).^2);
        funtI_Pt1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtI_St0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_St0(y1,y3)./(r(y1,y3).^2);
        funtI_St1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funI_St1(y1,y3)./(r(y1,y3).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_J_* 
        funtJ_Pt0_1_1 = @(y1,y3) (r1(y1,y3).^2).*funJ_Pt0(y1,y3)./(r(y1,y3).^2);
        funtJ_Pt1_1_1 = @(y1,y3) (r1(y1,y3).^2).*funJ_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtJ_St0_1_1 = @(y1,y3) (r1(y1,y3).^2).*funJ_St0(y1,y3)./(r(y1,y3).^2);
        funtJ_St1_1_1 = @(y1,y3) (r1(y1,y3).^2).*funJ_St1(y1,y3)./(r(y1,y3).^2);
        
        funtJ_Pt0_3_3 = @(y1,y3) (r3(y1,y3).^2).*funJ_Pt0(y1,y3)./(r(y1,y3).^2);
        funtJ_Pt1_3_3 = @(y1,y3) (r3(y1,y3).^2).*funJ_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtJ_St0_3_3 = @(y1,y3) (r3(y1,y3).^2).*funJ_St0(y1,y3)./(r(y1,y3).^2);
        funtJ_St1_3_3 = @(y1,y3) (r3(y1,y3).^2).*funJ_St1(y1,y3)./(r(y1,y3).^2);
     
        funtJ_Pt0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_Pt0(y1,y3)./(r(y1,y3).^2);
        funtJ_Pt1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_Pt1(y1,y3)./(r(y1,y3).^2);
        
        funtJ_St0_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_St0(y1,y3)./(r(y1,y3).^2);
        funtJ_St1_1_3 = @(y1,y3) (r1(y1,y3).*r3(y1,y3)).*funJ_St1(y1,y3)./(r(y1,y3).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni integrande   
        
        g_2 = @(y1,y3) +(funI_Pt0(y1,y3)-funI_St0(y1,y3))-funJ_St0(y1,y3)...
            -(funI_Pt1(y1,y3)-funI_St1(y1,y3))+funJ_St1(y1,y3);
        
        g_1 = @(y1,y3) -(funtJ_Pt0_1_1(y1,y3)-funtJ_St0_1_1(y1,y3))+(funI_Pt0(y1,y3)-funI_St0(y1,y3))-3*(funtI_Pt0_1_1(y1,y3)-funtI_St0_1_1(y1,y3))-funJ_St0(y1,y3)...
            +(funtJ_Pt1_1_1(y1,y3)-funtJ_St1_1_1(y1,y3))-(funI_Pt1(y1,y3)-funI_St1(y1,y3))+3*(funtI_Pt1_1_1(y1,y3)-funtI_St1_1_1(y1,y3))+funJ_St1(y1,y3)...
             -(funtJ_Pt0_1_3(y1,y3)-funtJ_St0_1_3(y1,y3))-3*(funtI_Pt0_1_3(y1,y3)-funtI_St0_1_3(y1,y3))...
            +(funtJ_Pt1_1_3(y1,y3)-funtJ_St1_1_3(y1,y3))+3*(funtI_Pt1_1_3(y1,y3)-funtI_St1_1_3(y1,y3));

        
        g_3 = @(y1,y3) -(funtJ_Pt0_1_3(y1,y3)-funtJ_St0_1_3(y1,y3))-3*(funtI_Pt0_1_3(y1,y3)-funtI_St0_1_3(y1,y3))...
            +(funtJ_Pt1_1_3(y1,y3)-funtJ_St1_1_3(y1,y3))+3*(funtI_Pt1_1_3(y1,y3)-funtI_St1_1_3(y1,y3))...
            -(funtJ_Pt0_3_3(y1,y3)-funtJ_St0_3_3(y1,y3))+(funI_Pt0(y1,y3)-funI_St0(y1,y3))-3*(funtI_Pt0_3_3(y1,y3)-funtI_St0_3_3(y1,y3))-funJ_St0(y1,y3)...
            +(funtJ_Pt1_3_3(y1,y3)-funtJ_St1_3_3(y1,y3))-(funI_Pt1(y1,y3)-funI_St1(y1,y3))+3*(funtI_Pt1_3_3(y1,y3)-funtI_St1_3_3(y1,y3))+funJ_St1(y1,y3);
                     
        ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        ris(2) = integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        ris(3) = integral2(g_3,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        
        ris(1) = ris(1)/(4*pi*pb_param.rho);
        ris(2) = ris(2)/(4*pi*pb_param.rho);
        ris(3) = ris(3)/(4*pi*pb_param.rho);
        
end
 
return

% function ris = time3D_intLE_rhs(L,pb_param,x1,x2,t0,t1)
% %-------------------------------------------------------------------------
% 
% %Inizializzazione dei vettori contenenti i RISULTATI dell'INTEGRAZIONE
% ris = zeros(3,1); 
% 
% velC_S = pb_param.velS;
% velC_P = pb_param.velP;
% 
% r = @(y1,y3) sqrt((x1-y1).^2+(x2-y3).^2);
% 
% fun1_Pt0 = @(y1,y3) heaviside(velC_P*t0-r(y1,y3)).*(t0-r(y1,y3)/velC_P).*(t0+r(y1,y3)/velC_P)./(2*r(y1,y3).^3);
% fun1_Pt1 = @(y1,y3) heaviside(velC_P*t1-r(y1,y3)).*(t1-r(y1,y3)/velC_P).*(t1+r(y1,y3)/velC_P)./(2*r(y1,y3).^3);
% 
% fun1_St0 = @(y1,y3) heaviside(velC_S*t0-r(y1,y3)).*(t0-r(y1,y3)/velC_S).*(t0+r(y1,y3)/velC_S)./(2*r(y1,y3).^3);
% fun1_St1 = @(y1,y3) heaviside(velC_S*t1-r(y1,y3)).*(t1-r(y1,y3)/velC_S).*(t1+r(y1,y3)/velC_S)./(2*r(y1,y3).^3);
% 
% fun2_St0 = @(y1,y3) -heaviside(velC_S*t0-r(y1,y3))./(velC_S^2*r(y1,y3));
% fun2_St1 = @(y1,y3) -heaviside(velC_S*t1-r(y1,y3))./(velC_S^2*r(y1,y3));
% 
% fun = @(y1,y3) (fun1_Pt0(y1,y3)-fun1_St0(y1,y3)+fun2_St0(y1,y3))-...
%                (fun1_Pt1(y1,y3)-fun1_St1(y1,y3)+fun2_St1(y1,y3));
% 
% ris(3) = integral2(fun,-L,L,-L,L);
% 
% %  y = @(y1) -y1+1;
% %  
% %  ris(3) = integral2(fun,0,1,y,1);
% 
% return

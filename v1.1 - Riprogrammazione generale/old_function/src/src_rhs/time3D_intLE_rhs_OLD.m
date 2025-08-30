function ris = time3D_intLE_rhs(L,pb_param,x1,x2,t0,t1,ind_RHS)

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
r1 = @(y1,y2) y1-x1;
r2 = @(y1,y2) y2-x2;

%Definiamo il MODULO del VETTORE DISTANZA r tra un PUNTO di CAMPO
%e il PUNTO SORGENTE fissato (nodo di Gauss-Hammer)
r = @(y1,y2) sqrt(r1(y1,y2).^2+r2(y1,y2).^2);

switch ind_RHS
    
    case 1 
        %Dato che permette di ottenere come soluzione phi = (1,0,0) sulla
        %screen quadrata
 
        %----------------------------------------------------------------
        %Definizione delle funzioni I_*
        funI_Pt0 = @(y1,y2) (heaviside(t0-r(y1,y2)/velC_P).*(t0-r(y1,y2)/velC_P).*(t0+r(y1,y2)/velC_P))./(2*r(y1,y2).^3);
        funI_Pt1 = @(y1,y2) (heaviside(t1-r(y1,y2)/velC_P).*(t1-r(y1,y2)/velC_P).*(t1+r(y1,y2)/velC_P))./(2*r(y1,y2).^3);
        
        funI_St0 = @(y1,y2) (heaviside(t0-r(y1,y2)/velC_S).*(t0-r(y1,y2)/velC_S).*(t0+r(y1,y2)/velC_S))./(2*r(y1,y2).^3);
        funI_St1 = @(y1,y2) (heaviside(t1-r(y1,y2)/velC_S).*(t1-r(y1,y2)/velC_S).*(t1+r(y1,y2)/velC_S))./(2*r(y1,y2).^3);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni J_*
        funJ_St0 = @(y1,y2) heaviside(t0-r(y1,y2)/velC_S)./(velC_S^2*r(y1,y2));
        funJ_St1 = @(y1,y2) heaviside(t1-r(y1,y2)/velC_S)./(velC_S^2*r(y1,y2));
        
        funJ_Pt0 = @(y1,y2) heaviside(t0-r(y1,y2)/velC_P)./(velC_P^2*r(y1,y2));
        funJ_Pt1 = @(y1,y2) heaviside(t1-r(y1,y2)/velC_P)./(velC_P^2*r(y1,y2));
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_I_*
        
        funtI_Pt0_1_1 = @(y1,y2) (r1(y1,y2).^2).*funI_Pt0(y1,y2)./(r(y1,y2).^2);
        funtI_Pt1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funI_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtI_St0_1_1 = @(y1,y2) (r1(y1,y2).^2).*funI_St0(y1,y2)./(r(y1,y2).^2);
        funtI_St1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funI_St1(y1,y2)./(r(y1,y2).^2);
        
        funtI_Pt0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_Pt0(y1,y2)./(r(y1,y2).^2);
        funtI_Pt1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtI_St0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_St0(y1,y2)./(r(y1,y2).^2);
        funtI_St1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_St1(y1,y2)./(r(y1,y2).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_J_*
        
        funtJ_Pt0_1_1 = @(y1,y2) (r1(y1,y2).^2).*funJ_Pt0(y1,y2)./(r(y1,y2).^2);
        funtJ_Pt1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funJ_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtJ_St0_1_1 = @(y1,y2) (r1(y1,y2).^2).*funJ_St0(y1,y2)./(r(y1,y2).^2);
        funtJ_St1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funJ_St1(y1,y2)./(r(y1,y2).^2);
        
        funtJ_Pt0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_Pt0(y1,y2)./(r(y1,y2).^2);
        funtJ_Pt1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtJ_St0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_St0(y1,y2)./(r(y1,y2).^2);
        funtJ_St1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_St1(y1,y2)./(r(y1,y2).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni integrande
        g_1 = @(y1,y2) -(funtJ_Pt0_1_1(y1,y2)-funtJ_St0_1_1(y1,y2))+(funI_Pt0(y1,y2)-funI_St0(y1,y2))-3*(funtI_Pt0_1_1(y1,y2)-funtI_St0_1_1(y1,y2))-funJ_St0(y1,y2)...
            +(funtJ_Pt1_1_1(y1,y2)-funtJ_St1_1_1(y1,y2))-(funI_Pt1(y1,y2)-funI_St1(y1,y2))+3*(funtI_Pt1_1_1(y1,y2)-funtI_St1_1_1(y1,y2))+funJ_St1(y1,y2);

        g_2 = @(y1,y2) -(funtJ_Pt0_1_2(y1,y2)-funtJ_St0_1_2(y1,y2))-3*(funtI_Pt0_1_2(y1,y2)-funtI_St0_1_2(y1,y2))...
            +(funtJ_Pt1_1_2(y1,y2)-funtJ_St1_1_2(y1,y2))+3*(funtI_Pt1_1_2(y1,y2)-funtI_St1_1_2(y1,y2));
        %----------------------------------------------------------------
      

        %Calcolo degli integrali
        ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        ris(2) = integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%         ris(1) = integral2(g_1,-L,L,-L,L);
%         ris(2) = integral2(g_2,-L,L,-L,L);
        %----------------------------------------------------------------
              
        
    case 2
        %Dato che permette di ottenere come soluzione phi = (0,1,0) sulla
        %screen quadrata
        
        %----------------------------------------------------------------
        %Definizione delle funzioni I_*
        funI_Pt0 = @(y1,y2) (heaviside(t0-r(y1,y2)/velC_P).*(t0-r(y1,y2)/velC_P).*(t0+r(y1,y2)/velC_P))./(2*r(y1,y2).^3);
        funI_Pt1 = @(y1,y2) (heaviside(t1-r(y1,y2)/velC_P).*(t1-r(y1,y2)/velC_P).*(t1+r(y1,y2)/velC_P))./(2*r(y1,y2).^3);
        
        funI_St0 = @(y1,y2) (heaviside(t0-r(y1,y2)/velC_S).*(t0-r(y1,y2)/velC_S).*(t0+r(y1,y2)/velC_S))./(2*r(y1,y2).^3);
        funI_St1 = @(y1,y2) (heaviside(t1-r(y1,y2)/velC_S).*(t1-r(y1,y2)/velC_S).*(t1+r(y1,y2)/velC_S))./(2*r(y1,y2).^3);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni J_*
        funJ_St0 = @(y1,y2) heaviside(t0-r(y1,y2)/velC_S)./(velC_S^2*r(y1,y2));
        funJ_St1 = @(y1,y2) heaviside(t1-r(y1,y2)/velC_S)./(velC_S^2*r(y1,y2));
        
        funJ_Pt0 = @(y1,y2) heaviside(t0-r(y1,y2)/velC_P)./(velC_P^2*r(y1,y2));
        funJ_Pt1 = @(y1,y2) heaviside(t1-r(y1,y2)/velC_P)./(velC_P^2*r(y1,y2));
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_I_*
        
        funtI_Pt0_2_2 = @(y1,y2) (r2(y1,y2).^2).*funI_Pt0(y1,y2)./(r(y1,y2).^2);
        funtI_Pt1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funI_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtI_St0_2_2 = @(y1,y2) (r2(y1,y2).^2).*funI_St0(y1,y2)./(r(y1,y2).^2);
        funtI_St1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funI_St1(y1,y2)./(r(y1,y2).^2);
        
        funtI_Pt0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_Pt0(y1,y2)./(r(y1,y2).^2);
        funtI_Pt1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtI_St0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_St0(y1,y2)./(r(y1,y2).^2);
        funtI_St1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_St1(y1,y2)./(r(y1,y2).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_J_*
        
        funtJ_Pt0_2_2 = @(y1,y2) (r2(y1,y2).^2).*funJ_Pt0(y1,y2)./(r(y1,y2).^2);
        funtJ_Pt1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funJ_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtJ_St0_2_2 = @(y1,y2) (r2(y1,y2).^2).*funJ_St0(y1,y2)./(r(y1,y2).^2);
        funtJ_St1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funJ_St1(y1,y2)./(r(y1,y2).^2);
        
        funtJ_Pt0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_Pt0(y1,y2)./(r(y1,y2).^2);
        funtJ_Pt1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtJ_St0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_St0(y1,y2)./(r(y1,y2).^2);
        funtJ_St1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_St1(y1,y2)./(r(y1,y2).^2);
        %----------------------------------------------------------------
  
        %Definizione delle funzioni integrande
        g_1 = @(y1,y2) -(funtJ_Pt0_1_2(y1,y2)-funtJ_St0_1_2(y1,y2))-3*(funtI_Pt0_1_2(y1,y2)-funtI_St0_1_2(y1,y2))...
            +(funtJ_Pt1_1_2(y1,y2)-funtJ_St1_1_2(y1,y2))+3*(funtI_Pt1_1_2(y1,y2)-funtI_St1_1_2(y1,y2));

        g_2 = @(y1,y2) -(funtJ_Pt0_2_2(y1,y2)-funtJ_St0_2_2(y1,y2))+(funI_Pt0(y1,y2)-funI_St0(y1,y2))-3*(funtI_Pt0_2_2(y1,y2)-funtI_St0_2_2(y1,y2))-funJ_St0(y1,y2)...
            +(funtJ_Pt1_2_2(y1,y2)-funtJ_St1_2_2(y1,y2))-(funI_Pt1(y1,y2)-funI_St1(y1,y2))+3*(funtI_Pt1_2_2(y1,y2)-funtI_St1_2_2(y1,y2))+funJ_St1(y1,y2);
        %----------------------------------------------------------------
        %Calcolo degli integrali    
        ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-7,'RelTol',1e-7);
        ris(2) = integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-5,'RelTol',1e-5);
        %----------------------------------------------------------------
         
    case 3 
        %Dato che permette di ottenere come soluzione phi = (0,0,1) sulla
        %screen quadrata
        
        %----------------------------------------------------------------
        %Definizione delle funzioni I_*
        funI_Pt0 = @(y1,y2) (heaviside(t0-r(y1,y2)/velC_P).*(t0-r(y1,y2)/velC_P).*(t0+r(y1,y2)/velC_P))./(2*r(y1,y2).^3);
        funI_Pt1 = @(y1,y2) (heaviside(t1-r(y1,y2)/velC_P).*(t1-r(y1,y2)/velC_P).*(t1+r(y1,y2)/velC_P))./(2*r(y1,y2).^3);
        
        funI_St0 = @(y1,y2) (heaviside(t0-r(y1,y2)/velC_S).*(t0-r(y1,y2)/velC_S).*(t0+r(y1,y2)/velC_S))./(2*r(y1,y2).^3);
        funI_St1 = @(y1,y2) (heaviside(t1-r(y1,y2)/velC_S).*(t1-r(y1,y2)/velC_S).*(t1+r(y1,y2)/velC_S))./(2*r(y1,y2).^3);
        %----------------------------------------------------------------
        
        %Definizione delle funzioni J_*
        funJ_St0 = @(y1,y2) heaviside(t0-r(y1,y2)/velC_S)./(velC_S^2*r(y1,y2));
        funJ_St1 = @(y1,y2) heaviside(t1-r(y1,y2)/velC_S)./(velC_S^2*r(y1,y2));
        %----------------------------------------------------------------
        
        %Definizione della funzione integranda
         g_3 = @(y1,y2) +(funI_Pt0(y1,y2)-funI_St0(y1,y2))-funJ_St0(y1,y2)...
            -(funI_Pt1(y1,y2)-funI_St1(y1,y2))+funJ_St1(y1,y2);
        %----------------------------------------------------------------
        
        %Calcolo dell'integrale
        ris(3) = integral2(g_3,-L,L,-L,L,'Method','iterated','AbsTol',1e-5,'RelTol',1e-5);
        
    case 4
        %Dato che permette di ottenere come soluzione phi = (1,1,1) sulla
        %screen quadrata
        
        %----------------------------------------------------------------
        %Definizione delle funzioni I_*
        funI_Pt0 = @(y1,y2) (heaviside(t0-r(y1,y2)/velC_P).*(t0-r(y1,y2)/velC_P).*(t0+r(y1,y2)/velC_P))./(2*r(y1,y2).^3);
        funI_Pt1 = @(y1,y2) (heaviside(t1-r(y1,y2)/velC_P).*(t1-r(y1,y2)/velC_P).*(t1+r(y1,y2)/velC_P))./(2*r(y1,y2).^3);
        
        funI_St0 = @(y1,y2) (heaviside(t0-r(y1,y2)/velC_S).*(t0-r(y1,y2)/velC_S).*(t0+r(y1,y2)/velC_S))./(2*r(y1,y2).^3);
        funI_St1 = @(y1,y2) (heaviside(t1-r(y1,y2)/velC_S).*(t1-r(y1,y2)/velC_S).*(t1+r(y1,y2)/velC_S))./(2*r(y1,y2).^3);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni J_*
        funJ_St0 = @(y1,y2) heaviside(t0-r(y1,y2)/velC_S)./(velC_S^2*r(y1,y2));
        funJ_St1 = @(y1,y2) heaviside(t1-r(y1,y2)/velC_S)./(velC_S^2*r(y1,y2));
        
        funJ_Pt0 = @(y1,y2) heaviside(t0-r(y1,y2)/velC_P)./(velC_P^2*r(y1,y2));
        funJ_Pt1 = @(y1,y2) heaviside(t1-r(y1,y2)/velC_P)./(velC_P^2*r(y1,y2));
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_I_*   
        funtI_Pt0_1_1 = @(y1,y2) (r1(y1,y2).^2).*funI_Pt0(y1,y2)./(r(y1,y2).^2);
        funtI_Pt1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funI_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtI_St0_1_1 = @(y1,y2) (r1(y1,y2).^2).*funI_St0(y1,y2)./(r(y1,y2).^2);
        funtI_St1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funI_St1(y1,y2)./(r(y1,y2).^2);
        
        funtI_Pt0_2_2 = @(y1,y2) (r2(y1,y2).^2).*funI_Pt0(y1,y2)./(r(y1,y2).^2);
        funtI_Pt1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funI_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtI_St0_2_2 = @(y1,y2) (r2(y1,y2).^2).*funI_St0(y1,y2)./(r(y1,y2).^2);
        funtI_St1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funI_St1(y1,y2)./(r(y1,y2).^2);
        
        funtI_Pt0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_Pt0(y1,y2)./(r(y1,y2).^2);
        funtI_Pt1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtI_St0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_St0(y1,y2)./(r(y1,y2).^2);
        funtI_St1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_St1(y1,y2)./(r(y1,y2).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni tilde_J_* 
        funtJ_Pt0_1_1 = @(y1,y2) (r1(y1,y2).^2).*funJ_Pt0(y1,y2)./(r(y1,y2).^2);
        funtJ_Pt1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funJ_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtJ_St0_1_1 = @(y1,y2) (r1(y1,y2).^2).*funJ_St0(y1,y2)./(r(y1,y2).^2);
        funtJ_St1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funJ_St1(y1,y2)./(r(y1,y2).^2);
        
        funtJ_Pt0_2_2 = @(y1,y2) (r2(y1,y2).^2).*funJ_Pt0(y1,y2)./(r(y1,y2).^2);
        funtJ_Pt1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funJ_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtJ_St0_2_2 = @(y1,y2) (r2(y1,y2).^2).*funJ_St0(y1,y2)./(r(y1,y2).^2);
        funtJ_St1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funJ_St1(y1,y2)./(r(y1,y2).^2);
     
        funtJ_Pt0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_Pt0(y1,y2)./(r(y1,y2).^2);
        funtJ_Pt1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_Pt1(y1,y2)./(r(y1,y2).^2);
        
        funtJ_St0_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_St0(y1,y2)./(r(y1,y2).^2);
        funtJ_St1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_St1(y1,y2)./(r(y1,y2).^2);
        %----------------------------------------------------------------
         
        %Definizione delle funzioni integrande   
        
        g_1 = @(y1,y2) -(funtJ_Pt0_1_1(y1,y2)-funtJ_St0_1_1(y1,y2))+(funI_Pt0(y1,y2)-funI_St0(y1,y2))-3*(funtI_Pt0_1_1(y1,y2)-funtI_St0_1_1(y1,y2))-funJ_St0(y1,y2)...
            +(funtJ_Pt1_1_1(y1,y2)-funtJ_St1_1_1(y1,y2))-(funI_Pt1(y1,y2)-funI_St1(y1,y2))+3*(funtI_Pt1_1_1(y1,y2)-funtI_St1_1_1(y1,y2))+funJ_St1(y1,y2)...
             -(funtJ_Pt0_1_2(y1,y2)-funtJ_St0_1_2(y1,y2))-3*(funtI_Pt0_1_2(y1,y2)-funtI_St0_1_2(y1,y2))...
            +(funtJ_Pt1_1_2(y1,y2)-funtJ_St1_1_2(y1,y2))+3*(funtI_Pt1_1_2(y1,y2)-funtI_St1_1_2(y1,y2));

        
        g_2 = @(y1,y2) -(funtJ_Pt0_1_2(y1,y2)-funtJ_St0_1_2(y1,y2))-3*(funtI_Pt0_1_2(y1,y2)-funtI_St0_1_2(y1,y2))...
            +(funtJ_Pt1_1_2(y1,y2)-funtJ_St1_1_2(y1,y2))+3*(funtI_Pt1_1_2(y1,y2)-funtI_St1_1_2(y1,y2))...
            -(funtJ_Pt0_2_2(y1,y2)-funtJ_St0_2_2(y1,y2))+(funI_Pt0(y1,y2)-funI_St0(y1,y2))-3*(funtI_Pt0_2_2(y1,y2)-funtI_St0_2_2(y1,y2))-funJ_St0(y1,y2)...
            +(funtJ_Pt1_2_2(y1,y2)-funtJ_St1_2_2(y1,y2))-(funI_Pt1(y1,y2)-funI_St1(y1,y2))+3*(funtI_Pt1_2_2(y1,y2)-funtI_St1_2_2(y1,y2))+funJ_St1(y1,y2);
                     
        g_3 = @(y1,y2) +(funI_Pt0(y1,y2)-funI_St0(y1,y2))-funJ_St0(y1,y2)...
            -(funI_Pt1(y1,y2)-funI_St1(y1,y2))+funJ_St1(y1,y2);
        
        ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-5,'RelTol',1e-5);
        ris(2) = integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-5,'RelTol',1e-5);
        ris(3) = integral2(g_3,-L,L,-L,L,'Method','iterated','AbsTol',1e-5,'RelTol',1e-5);
        
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
% r = @(y1,y2) sqrt((x1-y1).^2+(x2-y2).^2);
% 
% fun1_Pt0 = @(y1,y2) heaviside(velC_P*t0-r(y1,y2)).*(t0-r(y1,y2)/velC_P).*(t0+r(y1,y2)/velC_P)./(2*r(y1,y2).^3);
% fun1_Pt1 = @(y1,y2) heaviside(velC_P*t1-r(y1,y2)).*(t1-r(y1,y2)/velC_P).*(t1+r(y1,y2)/velC_P)./(2*r(y1,y2).^3);
% 
% fun1_St0 = @(y1,y2) heaviside(velC_S*t0-r(y1,y2)).*(t0-r(y1,y2)/velC_S).*(t0+r(y1,y2)/velC_S)./(2*r(y1,y2).^3);
% fun1_St1 = @(y1,y2) heaviside(velC_S*t1-r(y1,y2)).*(t1-r(y1,y2)/velC_S).*(t1+r(y1,y2)/velC_S)./(2*r(y1,y2).^3);
% 
% fun2_St0 = @(y1,y2) -heaviside(velC_S*t0-r(y1,y2))./(velC_S^2*r(y1,y2));
% fun2_St1 = @(y1,y2) -heaviside(velC_S*t1-r(y1,y2))./(velC_S^2*r(y1,y2));
% 
% fun = @(y1,y2) (fun1_Pt0(y1,y2)-fun1_St0(y1,y2)+fun2_St0(y1,y2))-...
%                (fun1_Pt1(y1,y2)-fun1_St1(y1,y2)+fun2_St1(y1,y2));
% 
% ris(3) = integral2(fun,-L,L,-L,L);
% 
% %  y = @(y1) -y1+1;
% %  
% %  ris(3) = integral2(fun,0,1,y,1);
% 
% return

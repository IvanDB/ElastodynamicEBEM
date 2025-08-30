function ris = time3D_intLE_rhs(L,pb_param,sp,t0,t1,ind_RHS)

%La function time3D_intLE() prende in INPUT:
%- la variabile L contenente la lunghezza del lato della screen quadrata
%
%- la struct pb_param contenente i parametri del problema
%
%- il vettore 3x1 contenente le coordinate del punto sorgente sp
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

%Ricaviamo l'ASCISSA e l'ORDINATA del PUNTO SORGENTE
x1=sp(1);
x2=sp(2);

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
 
%************************************************************************               
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
%         %----------------------------------------------------------------
%          
        %Definizione delle funzioni integrande
        g_1 = @(y1,y2) -(funtJ_Pt0_1_1(y1,y2)-funtJ_St0_1_1(y1,y2))+(funI_Pt0(y1,y2)-funI_St0(y1,y2))-3*(funtI_Pt0_1_1(y1,y2)-funtI_St0_1_1(y1,y2))-funJ_St0(y1,y2)...
            +(funtJ_Pt1_1_1(y1,y2)-funtJ_St1_1_1(y1,y2))-(funI_Pt1(y1,y2)-funI_St1(y1,y2))+3*(funtI_Pt1_1_1(y1,y2)-funtI_St1_1_1(y1,y2))+funJ_St1(y1,y2);

        g_2 = @(y1,y2) -(funtJ_Pt0_1_2(y1,y2)-funtJ_St0_1_2(y1,y2))-3*(funtI_Pt0_1_2(y1,y2)-funtI_St0_1_2(y1,y2))...
             +(funtJ_Pt1_1_2(y1,y2)-funtJ_St1_1_2(y1,y2))+3*(funtI_Pt1_1_2(y1,y2)-funtI_St1_1_2(y1,y2));
        %Calcolo degli integrali
        ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        ris(2) = integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%************************************************************************       

%         g_1 = @(y1,y2) -(funtJ_Pt0_1_1(y1,y2)-funtJ_St0_1_1(y1,y2))+(funtJ_Pt1_1_1(y1,y2)-funtJ_St1_1_1(y1,y2));
%         ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
%         
%         g_1 = @(y1,y2)-funJ_St0(y1,y2)+funJ_St1(y1,y2);
%         ris(1) = ris(1)+ integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
%         
%         g_1 = @(y1,y2) +(funI_Pt0(y1,y2)-funI_St0(y1,y2))-(funI_Pt1(y1,y2)-funI_St1(y1,y2));
%         ris(1) = ris(1)+ integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
% 
%         g_1 = @(y1,y2) -3*(funtI_Pt0_1_1(y1,y2)-funtI_St0_1_1(y1,y2))+3*(funtI_Pt1_1_1(y1,y2)-funtI_St1_1_1(y1,y2));
%         ris(1) = ris(1)+ integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
% 
%         g_2 = @(y1,y2) -(funtJ_Pt0_1_2(y1,y2)-funtJ_St0_1_2(y1,y2))...
%                         +(funtJ_Pt1_1_2(y1,y2)-funtJ_St1_1_2(y1,y2));
%         ris(2) = integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
% 
%         g_2 = @(y1,y2)-3*(funtI_Pt0_1_2(y1,y2)-funtI_St0_1_2(y1,y2))+3*(funtI_Pt1_1_2(y1,y2)-funtI_St1_1_2(y1,y2));
%         ris(2) = ris(2)+integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);    
%        ----------------------------------------------------------------
%************************************************************************               
        
%           C1=(velC_P^2+velC_S^2)/(2*velC_P^2*velC_S^2);
%           C2=(velC_P^2-velC_S^2)/(2*velC_P^2*velC_S^2);
%           
%           H0_1=@(y1,y2) and((r(y1,y2)<=velC_S*t0), (r(y1,y2)<=velC_P*t0));
%           H1_1=@(y1,y2) and((r(y1,y2)<=velC_S*t1), (r(y1,y2)<=velC_P*t1));    
%           H0_2=@(y1,y2) and((r(y1,y2)>velC_S*t0),(r(y1,y2)<=velC_P*t0));
%           H1_2=@(y1,y2) and((r(y1,y2)>velC_S*t1), (r(y1,y2)<=velC_P*t1));
%           
%           %Primo pezzo
%           g_1_1=@(y1,y2) C1./r(y1,y2);
%           gt_1_1=@(y1,y2) (-H0_1(y1,y2)+H1_1(y1,y2)).*g_1_1(y1,y2);
%           ris(1) = integral2(gt_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%           
%           g_1_1=@(y1,y2) +C2*((r1(y1,y2).^2)./(r(y1,y2).^3));
%           gt_1_1=@(y1,y2) (-H0_1(y1,y2)+H1_1(y1,y2)).*g_1_1(y1,y2);
%           ris(1) = ris(1)+integral2(gt_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%    
%           %Secondo pezzo
%           g_1_1=@(y1,y2) 1./(2*velC_P^2*r(y1,y2));
%           gt_1_1=@(y1,y2) (-H0_2(y1,y2)+H1_2(y1,y2)).*g_1_1(y1,y2);
%           ris(1) = ris(1)+ integral2(gt_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%           g_1_1=@(y1,y2) -(r1(y1,y2).^2)./(2*velC_P^2*r(y1,y2).^3);
%           gt_1_1=@(y1,y2) (-H0_2(y1,y2)+H1_2(y1,y2)).*g_1_1(y1,y2);
%           ris(1) = ris(1)+ integral2(gt_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%        
%           g_1_1=@(y1,y2) -1./(2*r(y1,y2).^3);
%           gt_1_1=@(y1,y2) (-t0^2*H0_2(y1,y2)+t1^2*H1_2(y1,y2)).*g_1_1(y1,y2);
%           ris(1) = ris(1)+ integral2(gt_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%           g_1_1=@(y1,y2) ((3*r1(y1,y2).^2)./(2*r(y1,y2).^5));
%           gt_1_1=@(y1,y2) (-t0^2*H0_2(y1,y2)+t1^2*H1_2(y1,y2)).*g_1_1(y1,y2);
%           ris(1) = ris(1)+ integral2(gt_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
% 
%           %Primo pezzo        
%           g_1_2=@(y1,y2) C2*((r1(y1,y2).*r2(y1,y2))./(r(y1,y2).^3));
%           gt_1_2=@(y1,y2) (-H0_1(y1,y2)+H1_1(y1,y2)).*g_1_2(y1,y2);
%           ris(2) = integral2(gt_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%           %Secondo pezzo
%           g_1_2=@(y1,y2) -(r1(y1,y2).*r2(y1,y2))./(2*velC_P^2*r(y1,y2).^3);
%           gt_1_2=@(y1,y2) (-H0_2(y1,y2)+H1_2(y1,y2)).*g_1_2(y1,y2);
%           ris(2) = ris(2) + integral2(gt_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%       
%           g_1_2=@(y1,y2) ((3*r1(y1,y2).*r2(y1,y2))./(2*r(y1,y2).^5));
%           gt_1_2=@(y1,y2) (-t0^2*H0_2(y1,y2)+t1^2*H1_2(y1,y2)).*g_1_2(y1,y2);
%           ris(2) = ris(2)+ integral2(gt_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%  
          %Dividiamo i risultati ottenuti per il coefficiente davnti al
          %nucleo
          ris(1) = ris(1)/(4*pi*pb_param.rho);
          ris(2) = ris(2)/(4*pi*pb_param.rho);

        %----------------------------------------------------------------
              
        
    case 2
        %Dato che permette di ottenere come soluzione phi = (0,1,0) sulla
        %screen quadrata
        
%************************************************************************              
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
  
      g_1 = @(y1,y2) -(funtJ_Pt0_1_2(y1,y2)-funtJ_St0_1_2(y1,y2))-3*(funtI_Pt0_1_2(y1,y2)-funtI_St0_1_2(y1,y2))...
            +(funtJ_Pt1_1_2(y1,y2)-funtJ_St1_1_2(y1,y2))+3*(funtI_Pt1_1_2(y1,y2)-funtI_St1_1_2(y1,y2));
      g_2 = @(y1,y2) -(funtJ_Pt0_2_2(y1,y2)-funtJ_St0_2_2(y1,y2))+(funI_Pt0(y1,y2)-funI_St0(y1,y2))-3*(funtI_Pt0_2_2(y1,y2)-funtI_St0_2_2(y1,y2))-funJ_St0(y1,y2)...
            +(funtJ_Pt1_2_2(y1,y2)-funtJ_St1_2_2(y1,y2))-(funI_Pt1(y1,y2)-funI_St1(y1,y2))+3*(funtI_Pt1_2_2(y1,y2)-funtI_St1_2_2(y1,y2))+funJ_St1(y1,y2);

% %-----------------------------------------------------------------------
        %Calcolo degli integrali    
        ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
        ris(2) = integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
% %----------------------------------------------------------------
%************************************************************************       


%Definizione delle funzioni integrande e calcolo degli integrali

%         g_1 = @(y1,y2) -(funtJ_Pt0_1_2(y1,y2)-funtJ_St0_1_2(y1,y2))+(funtJ_Pt1_1_2(y1,y2)-funtJ_St1_1_2(y1,y2));
%         ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
%        
%         g_1 = @(y1,y2)-3*(funtI_Pt0_1_2(y1,y2)-funtI_St0_1_2(y1,y2)) +3*(funtI_Pt1_1_2(y1,y2)-funtI_St1_1_2(y1,y2));
%         ris(1) = ris(1)+ integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
% 
%           g_2 = @(y1,y2) -(funtJ_Pt0_2_2(y1,y2)-funtJ_St0_2_2(y1,y2))+(funtJ_Pt1_2_2(y1,y2)-funtJ_St1_2_2(y1,y2));  
%           ris(2) = integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
% 
%           g_2 = @(y1,y2) -funJ_St0(y1,y2)+funJ_St1(y1,y2);
%           ris(2) = ris(2)+ integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
% 
%           g_2 = @(y1,y2) +(funI_Pt0(y1,y2)-funI_St0(y1,y2))-(funI_Pt1(y1,y2)-funI_St1(y1,y2));       
%           ris(2) = ris(2)+ integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
%             
%           g_2 = @(y1,y2)-3*(funtI_Pt0_2_2(y1,y2)-funtI_St0_2_2(y1,y2))+3*(funtI_Pt1_2_2(y1,y2)-funtI_St1_2_2(y1,y2));
%           ris(2) = ris(2)+ integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
%************************************************************************       


%           C1=(velC_P^2+velC_S^2)/(2*velC_P^2*velC_S^2);
%           C2=(velC_P^2-velC_S^2)/(2*velC_P^2*velC_S^2);
%           
%           H0_1=@(y1,y2) and((r(y1,y2)<=velC_S*t0), (r(y1,y2)<=velC_P*t0));
%           H1_1=@(y1,y2) and((r(y1,y2)<=velC_S*t1), (r(y1,y2)<=velC_P*t1));    
%           H0_2=@(y1,y2) and((r(y1,y2)>velC_S*t0),(r(y1,y2)<=velC_P*t0));
%           H1_2=@(y1,y2) and((r(y1,y2)>velC_S*t1), (r(y1,y2)<=velC_P*t1));
%           
%           %Primo pezzo        
%           g_1_2=@(y1,y2) C2*((r1(y1,y2).*r2(y1,y2))./(r(y1,y2).^3));
%           gt_1_2=@(y1,y2) (-H0_1(y1,y2)+H1_1(y1,y2)).*g_1_2(y1,y2);
%           ris(1) = integral2(gt_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%           %Secondo pezzo
%           g_1_2=@(y1,y2) -(r1(y1,y2).*r2(y1,y2))./(2*velC_P^2*r(y1,y2).^3);
%           gt_1_2=@(y1,y2) (-H0_2(y1,y2)+H1_2(y1,y2)).*g_1_2(y1,y2);
%           ris(1) = ris(1) + integral2(gt_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%       
%           g_1_2=@(y1,y2,t) ((3*r1(y1,y2).*r2(y1,y2))./(2*r(y1,y2).^5));
%           gt_1_2=@(y1,y2) (-t0^2.*H0_2(y1,y2)+t1^2.*H1_2(y1,y2)).*g_1_2(y1,y2,t0);
%           ris(1) = ris(1)+ integral2(gt_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%  
%           %Primo pezzo
%           g_2_2=@(y1,y2) C1./r(y1,y2);
%           gt_2_2=@(y1,y2) (-H0_1(y1,y2)+H1_1(y1,y2)).*g_2_2(y1,y2);
%           ris(2) = integral2(gt_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%           
%           g_2_2=@(y1,y2) +C2*((r2(y1,y2).^2)./(r(y1,y2).^3));
%           gt_2_2=@(y1,y2) (-H0_1(y1,y2)+H1_1(y1,y2)).*g_2_2(y1,y2);
%           ris(2) = ris(2)+integral2(gt_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%           
%           %Secondo pezzo
%           g_2_2=@(y1,y2) 1./(2*velC_P^2*r(y1,y2));
%           gt_2_2=@(y1,y2) (-H0_2(y1,y2)+H1_2(y1,y2)).*g_2_2(y1,y2);
%           ris(2) = ris(2)+ integral2(gt_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%           g_2_2=@(y1,y2) -(r2(y1,y2).^2)./(2*velC_P^2*r(y1,y2).^3);
%           gt_2_2=@(y1,y2) (-H0_2(y1,y2)+H1_2(y1,y2)).*g_2_2(y1,y2);
%           ris(2) = ris(2)+ integral2(gt_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%        
%           g_2_2=@(y1,y2) -1./(2*r(y1,y2).^3);
%           gt_2_2=@(y1,y2) (-t0^2*H0_2(y1,y2)+t1^2*H1_2(y1,y2)).*g_2_2(y1,y2);
%           ris(2) = ris(2)+ integral2(gt_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%           g_2_2=@(y1,y2) (3*r2(y1,y2).^2)./(2*r(y1,y2).^5);
%           gt_2_2=@(y1,y2) (-t0^2*H0_2(y1,y2)+t1^2*H1_2(y1,y2)).*g_2_2(y1,y2);
%           ris(2) = ris(2)+ integral2(gt_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);

          %Dividiamo i risultati ottenuti per il coefficiente davnti al
          %nucleo
          ris(1) = ris(1)/(4*pi*pb_param.rho);
          ris(2) = ris(2)/(4*pi*pb_param.rho);
         
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
        ris(3) = integral2(g_3,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        
        ris(3) = ris(3)/(4*pi*pb_param.rho);
        
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
        
        ris(1) = integral2(g_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        ris(2) = integral2(g_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        ris(3) = integral2(g_3,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-5);
        
        ris(1) = ris(1)/(4*pi*pb_param.rho);
        ris(2) = ris(2)/(4*pi*pb_param.rho);
        ris(3) = ris(3)/(4*pi*pb_param.rho);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     warning('off','all')
%     
%     C1=(velC_P^2+velC_S^2)/(2*velC_P^2*velC_S^2);
%     C2=(velC_P^2-velC_S^2)/(2*velC_P^2*velC_S^2);
% 
%     H0_1=@(y1,y2) and((r(y1,y2)<=velC_S*t0), (r(y1,y2)<=velC_P*t0));
%     H1_1=@(y1,y2) and((r(y1,y2)<=velC_S*t1), (r(y1,y2)<=velC_P*t1));    
%     H0_2=@(y1,y2) and((r(y1,y2)>velC_S*t0),(r(y1,y2)<=velC_P*t0));
%     H1_2=@(y1,y2) and((r(y1,y2)>velC_S*t1), (r(y1,y2)<=velC_P*t1));
%     
%     %Integrali utili
%     %---------------------------------------------------------------------
%     f_1_r=@(y1,y2) (-H0_1(y1,y2)+H1_1(y1,y2))./r(y1,y2);
%     int_1_r=integral2(f_1_r,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_1_r=int_1_r+integral2(f_1_r,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_1_r=int_1_r+integral2(f_1_r,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_1_r=int_1_r+integral2(f_1_r,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
% 
%     f2_1_r=@(y1,y2) (-H0_2(y1,y2)+H1_2(y1,y2))./r(y1,y2);
%     int2_1_r=integral2(f2_1_r,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int2_1_r=int2_1_r+integral2(f2_1_r,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int2_1_r=int2_1_r+integral2(f2_1_r,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int2_1_r=int2_1_r+integral2(f2_1_r,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
% 
%     f_1_r_3=@(y1,y2) (-t0^2*H0_2(y1,y2)+t1^2*H1_2(y1,y2))./(r(y1,y2).^3);
%     int_1_r_3=integral2(f_1_r_3,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_1_r_3=int_1_r_3+integral2(f_1_r_3,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_1_r_3=int_1_r_3+integral2(f_1_r_3,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_1_r_3=int_1_r_3+integral2(f_1_r_3,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     %---------------------------------------------------------------------
%     
%     %---------------------------------------------------------------------
%     %INTEGRALE di g_1_2=g_2_1
%     %---------------------------------------------------------------------
%     %Primo pezzo g_1_2=g_2_1
%     g_1_2=@(y1,y2) C2*(-H0_1(y1,y2)+H1_1(y1,y2)).*((r1(y1,y2).*r2(y1,y2))./(r(y1,y2).^3));
%     int_g_1_2=integral2(g_1_2,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_g_1_2=int_g_1_2+integral2(g_1_2,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_g_1_2=int_g_1_2+integral2(g_1_2,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_g_1_2=int_g_1_2+integral2(g_1_2,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
% 
%     
%     %Secondo pezzo di g_1_2=g_2_1
%     g_1_2=@(y1,y2) (-H0_2(y1,y2)+H1_2(y1,y2)).*(-(r1(y1,y2).*r2(y1,y2))./(2*velC_P^2*r(y1,y2).^3));
%     int_g_1_2=int_g_1_2+integral2(g_1_2,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_g_1_2=int_g_1_2+integral2(g_1_2,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_g_1_2=int_g_1_2+integral2(g_1_2,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_g_1_2=int_g_1_2+integral2(g_1_2,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
% 
%     g_1_2=@(y1,y2) (-t0^2*H0_2(y1,y2)+t1^2*H1_2(y1,y2)).*((3*r1(y1,y2).*r2(y1,y2))./(2*r(y1,y2).^5));
%     int_g_1_2=int_g_1_2+integral2(g_1_2,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_g_1_2=int_g_1_2+integral2(g_1_2,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_g_1_2=int_g_1_2+integral2(g_1_2,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     int_g_1_2=int_g_1_2+integral2(g_1_2,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     %---------------------------------------------------------------------
% 
%     %---------------------------------------------------------------------
%     %COSTRUZIONE DEL PRIMO ELEMENTO DEL TERMINE NOTO 
%     %---------------------------------------------------------------------
%     %Primo pezzo di g_1_1
%     ris(1) = C1*int_1_r;
% 
%     g_1_1=@(y1,y2) C2*(-H0_1(y1,y2)+H1_1(y1,y2)).*((r1(y1,y2).^2)./(r(y1,y2).^3));
%     ris(1)=ris(1)+integral2(g_1_1,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(1)=ris(1)+integral2(g_1_1,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(1)=ris(1)+integral2(g_1_1,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(1)=ris(1)+integral2(g_1_1,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
% 
%     %Secondo pezzo di g_1_1
%     ris(1) = ris(1)+ 1/(2*velC_P^2)*int2_1_r;
% 
%     g_1_1=@(y1,y2) (-H0_2(y1,y2)+H1_2(y1,y2)).*(-(r1(y1,y2).^2)./(2*velC_P^2*r(y1,y2).^3));
%     ris(1)=ris(1)+integral2(g_1_1,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(1)=ris(1)+integral2(g_1_1,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(1)=ris(1)+integral2(g_1_1,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(1)=ris(1)+integral2(g_1_1,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
% 
%     ris(1) = ris(1) -1/2*int_1_r_3;
% 
%     g_1_1=@(y1,y2) (-t0^2*H0_2(y1,y2)+t1^2*H1_2(y1,y2)).*((3*r1(y1,y2).^2)./(2*r(y1,y2).^5));
%     ris(1)=ris(1)+integral2(g_1_1,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(1)=ris(1)+integral2(g_1_1,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(1)=ris(1)+integral2(g_1_1,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(1)=ris(1)+integral2(g_1_1,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     
%     ris(1)=ris(1)+int_g_1_2;
%     %---------------------------------------------------------------------
%     
%     %---------------------------------------------------------------------
%     %COSTRUZIONE DEL SECONDO ELEMENTO DEL TERMINE NOTO 
%     %---------------------------------------------------------------------
%     %Primo pezzo di g_2_2
%     ris(2) = C1*int_1_r;
% 
%     g_2_2=@(y1,y2) C2*(-H0_1(y1,y2)+H1_1(y1,y2)).*((r2(y1,y2).^2)./(r(y1,y2).^3));
%     ris(2)=ris(2)+integral2(g_2_2,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(2)=ris(2)+integral2(g_2_2,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(2)=ris(2)+integral2(g_2_2,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(2)=ris(2)+integral2(g_2_2,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
% 
%     %Secondo pezzo di g_2_2
%     ris(2) = ris(2)+ 1./(2*velC_P^2)*int2_1_r;
% 
%     g_2_2=@(y1,y2) (-H0_2(y1,y2)+H1_2(y1,y2)).*(-(r2(y1,y2).^2)./(2*velC_P^2*r(y1,y2).^3));
%     ris(2)=ris(2)+integral2(g_2_2,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(2)=ris(2)+integral2(g_2_2,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(2)=ris(2)+integral2(g_2_2,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(2)=ris(2)+integral2(g_2_2,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
% % 
%     ris(2) = ris(2)-1/2*int_1_r_3;
% 
%     g_2_2=@(y1,y2) (-t0^2*H0_2(y1,y2)+t1^2*H1_2(y1,y2)).*((3*r2(y1,y2).^2)./(2*r(y1,y2).^5));
%     ris(2)=ris(2)+integral2(g_2_2,-L,x1,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(2)=ris(2)+integral2(g_2_2,x1,L,-L,x2,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(2)=ris(2)+integral2(g_2_2,-L,x1,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
%     ris(2)=ris(2)+integral2(g_2_2,x1,L,x2,L,'Method','iterated','RelTol',1e-6,'AbsTol',1e-6);
% 
%     ris(2)=ris(2)+int_g_1_2;
%     %---------------------------------------------------------------------
% 
%     %---------------------------------------------------------------------
%     %COSTRUZIONE DEL TERZO ELEMENTO DEL TERMINE NOTO 
%     %---------------------------------------------------------------------
%     %Primo pezzo di g_3_3
%     ris(3) =  C1*int_1_r;
% 
%     %Secondo pezzo di g_3_3
%     ris(3) = ris(3)+ 1./(2*velC_P^2)*int2_1_r;
% 
%     ris(3) = ris(3)-1/2*int_1_r_3;
% 
%     %---------------------------------------------------------------------
% 
%     ris(1) = ris(1)/(4*pi*pb_param.rho);
%     ris(2) = ris(2)/(4*pi*pb_param.rho);
%     ris(3) = ris(3)/(4*pi*pb_param.rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     warning('off','all')
%     
%     C1=(velC_P^2+velC_S^2)/(2*velC_P^2*velC_S^2);
%     C2=(velC_P^2-velC_S^2)/(2*velC_P^2*velC_S^2);
% 
%     H0_1=@(y1,y2) and((r(y1,y2)<=velC_S*t0), (r(y1,y2)<=velC_P*t0));
%     H1_1=@(y1,y2) and((r(y1,y2)<=velC_S*t1), (r(y1,y2)<=velC_P*t1));    
%     H0_2=@(y1,y2) and((r(y1,y2)>velC_S*t0),(r(y1,y2)<=velC_P*t0));
%     H1_2=@(y1,y2) and((r(y1,y2)>velC_S*t1), (r(y1,y2)<=velC_P*t1));
%     
%     %Integrali utili
%     %---------------------------------------------------------------------
%     f_1_r=@(y1,y2) (-H0_1(y1,y2))./r(y1,y2);
%     int_1_r=integral2(f_1_r,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%     f_1_r=@(y1,y2) (H1_1(y1,y2))./r(y1,y2);
%     int_1_r=int_1_r+integral2(f_1_r,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     
%     f2_1_r=@(y1,y2) (-H0_2(y1,y2))./r(y1,y2);
%     int2_1_r=integral2(f2_1_r,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     f2_1_r=@(y1,y2) (H1_2(y1,y2))./r(y1,y2);
%     int2_1_r=int2_1_r+integral2(f2_1_r,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%     f_1_r_3=@(y1,y2) (-t0^2*H0_2(y1,y2)+t1^2*H1_2(y1,y2))./(r(y1,y2).^3);
%     int_1_r_3=integral2(f_1_r_3,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% %     int_1_r_3=int_1_r_3+integral2(f_1_r_3,x1,L,-L,x2,'Method','tiled');
% %     int_1_r_3=int_1_r_3+integral2(f_1_r_3,-L,x1,x2,L,'Method','tiled');
% %     int_1_r_3=int_1_r_3+integral2(f_1_r_3,x1,L,x2,L,'Method','tiled');
%     %---------------------------------------------------------------------
%     
%     %---------------------------------------------------------------------
%     %INTEGRALE di g_1_2=g_2_1
%     %---------------------------------------------------------------------
%     %Primo pezzo g_1_2=g_2_1
%     g_1_2=@(y1,y2) C2*(-H0_1(y1,y2)).*((r1(y1,y2).*r2(y1,y2))./(r(y1,y2).^3));
%     int_g_1_2=integral2(g_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     g_1_2=@(y1,y2) C2*(+H1_1(y1,y2)).*((r1(y1,y2).*r2(y1,y2))./(r(y1,y2).^3));
%     int_g_1_2=int_g_1_2+integral2(g_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
% 
%     
%     %Secondo pezzo di g_1_2=g_2_1
%     g_1_2=@(y1,y2) (-H0_2(y1,y2)).*(-(r1(y1,y2).*r2(y1,y2))./(2*velC_P^2*r(y1,y2).^3));
%     int_g_1_2=int_g_1_2+integral2(g_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     g_1_2=@(y1,y2) (+H1_2(y1,y2)).*(-(r1(y1,y2).*r2(y1,y2))./(2*velC_P^2*r(y1,y2).^3));
%     int_g_1_2=int_g_1_2+integral2(g_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
% 
%     g_1_2=@(y1,y2) (-t0^2*H0_2(y1,y2)).*((3*r1(y1,y2).*r2(y1,y2))./(2*r(y1,y2).^5));
%     int_g_1_2=int_g_1_2+integral2(g_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     g_1_2=@(y1,y2) (+t1^2*H1_2(y1,y2)).*((3*r1(y1,y2).*r2(y1,y2))./(2*r(y1,y2).^5));
%     int_g_1_2=int_g_1_2+integral2(g_1_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%     %---------------------------------------------------------------------
% 
%     %---------------------------------------------------------------------
%     %COSTRUZIONE DEL PRIMO ELEMENTO DEL TERMINE NOTO 
%     %---------------------------------------------------------------------
%     %Primo pezzo di g_1_1
%     ris(1) = C1*int_1_r;
% 
%     g_1_1=@(y1,y2) C2*(-H0_1(y1,y2)).*((r1(y1,y2).^2)./(r(y1,y2).^3));
%     ris(1)=ris(1)+integral2(g_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     g_1_1=@(y1,y2) C2*(+H1_1(y1,y2)).*((r1(y1,y2).^2)./(r(y1,y2).^3));
%     ris(1)=ris(1)+integral2(g_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
% 
% %     g_1_1=@(y1,y2) -C2*H0_1(y1,y2).*((r1(y1,y2).^2)./(r(y1,y2).^3));
% %     ris(1)=ris(1)+integral2(g_1_1,-L,L,-L,L,'Method','tiled','AbsTol',1e-6,'RelTol',1e-6);
% % 
% %     g_1_1=@(y1,y2) C2*H1_1(y1,y2).*((r1(y1,y2).^2)./(r(y1,y2).^3));
% %     ris(1)=ris(1)+integral2(g_1_1,-L,L,-L,L,'Method','tiled','AbsTol',1e-6,'RelTol',1e-6);
% 
% 
%     %Secondo pezzo di g_1_1
%     ris(1) = ris(1)+ 1/(2*velC_P^2)*int2_1_r;
% 
%     g_1_1=@(y1,y2) (-H0_2(y1,y2)).*(-(r1(y1,y2).^2)./(2*velC_P^2*r(y1,y2).^3));
%     ris(1)=ris(1)+integral2(g_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     g_1_1=@(y1,y2) (+H1_2(y1,y2)).*(-(r1(y1,y2).^2)./(2*velC_P^2*r(y1,y2).^3));
%     ris(1)=ris(1)+integral2(g_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
% 
%     ris(1) = ris(1) -1/2*int_1_r_3;
% 
%     g_1_1=@(y1,y2) (-t0^2*H0_2(y1,y2)).*((3*r1(y1,y2).^2)./(2*r(y1,y2).^5));
%     ris(1)=ris(1)+integral2(g_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     g_1_1=@(y1,y2) (t1^2*H1_2(y1,y2)).*((3*r1(y1,y2).^2)./(2*r(y1,y2).^5));
%     ris(1)=ris(1)+integral2(g_1_1,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%     
%     ris(1)=ris(1)+int_g_1_2;
%     %---------------------------------------------------------------------
%     
%     %---------------------------------------------------------------------
%     %COSTRUZIONE DEL SECONDO ELEMENTO DEL TERMINE NOTO 
%     %---------------------------------------------------------------------
%     %Primo pezzo di g_2_2
%     ris(2) = C1*int_1_r;
% 
%     g_2_2=@(y1,y2) C2*(-H0_1(y1,y2)).*((r2(y1,y2).^2)./(r(y1,y2).^3));
%     ris(2)=ris(2)+integral2(g_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     g_2_2=@(y1,y2) C2*(+H1_1(y1,y2)).*((r2(y1,y2).^2)./(r(y1,y2).^3));
%     ris(2)=ris(2)+integral2(g_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%     %Secondo pezzo di g_2_2
%     ris(2) = ris(2)+ 1./(2*velC_P^2)*int2_1_r;
% 
%     g_2_2=@(y1,y2) (-H0_2(y1,y2)).*(-(r2(y1,y2).^2)./(2*velC_P^2*r(y1,y2).^3));
%     ris(2)=ris(2)+integral2(g_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     g_2_2=@(y1,y2) (+H1_2(y1,y2)).*(-(r2(y1,y2).^2)./(2*velC_P^2*r(y1,y2).^3));
%     ris(2)=ris(2)+integral2(g_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
% 
%     ris(2) = ris(2)-1/2*int_1_r_3;
% 
%     g_2_2=@(y1,y2) (-t0^2*H0_2(y1,y2)).*((3*r2(y1,y2).^2)./(2*r(y1,y2).^5));
%     ris(2)=ris(2)+integral2(g_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
%     g_2_2=@(y1,y2) (t1^2*H1_2(y1,y2)).*((3*r2(y1,y2).^2)./(2*r(y1,y2).^5));
%     ris(2)=ris(2)+integral2(g_2_2,-L,L,-L,L,'Method','iterated','AbsTol',1e-6,'RelTol',1e-6);
% 
%     ris(2)=ris(2)+int_g_1_2;
%     %---------------------------------------------------------------------
% 
%     %---------------------------------------------------------------------
%     %COSTRUZIONE DEL TERZO ELEMENTO DEL TERMINE NOTO 
%     %---------------------------------------------------------------------
%     %Primo pezzo di g_3_3
%     ris(3) =  C1*int_1_r;
% 
%     %Secondo pezzo di g_3_3
%     ris(3) = ris(3)+ 1./(2*velC_P^2)*int2_1_r;
% 
%     ris(3) = ris(3)-1/2*int_1_r_3;
% 
%     %---------------------------------------------------------------------
% 
%     ris(1) = ris(1)/(4*pi*pb_param.rho);
%     ris(2) = ris(2)/(4*pi*pb_param.rho);
%     ris(3) = ris(3)/(4*pi*pb_param.rho);
    
    case 5 
        
        if t0 <= 1.25e-01  
            %Se t0<=1/8
            
            ris(3) = ris(3) - ((sin(4*pi*t0)).^2)*x1;
            
        else
            %Se t0>1/8
            
            ris(3) = ris(3) - 1*x1;
            
        end
        
        if t1 <= 1.25e-01
            %Se t1<=1/8
            ris(3) = ris(3) + ((sin(4*pi*t1)).^2)*x1;
            
        else
            %Se t1>1/8
            ris(3) = ris(3) + 1*x1;
            
        end   
        
    case 6 %test bicchierino (1° versione)
        
        ris(3) = ((sin(t1))^5-(sin(t0))^5)*x1^2;
        
    case 7 %test Sfera
        
        ris(3) = (t1-t0)/2;
        
    case 8 %test barretts
        ris=test_bar(pb_param,sp,t0,t1);
         
        
    case 9 %Test bicchierino (2°versione)
        
        if t0 <= pi/2  
            %Se t0<=pi/2
            
            ris(3) = ris(3) - ((sin(t0))^5)*x1^2;
            
        else
            %Se t0>pi/2
            
            ris(3) = ris(3) - x1^2;
            
        end
        
        if t1 <= pi/2
            %Se t1<=pi/2
            ris(3) = ris(3) + ((sin(t1))^5)*x1^2;
            
        else
            %Se t1>pi/2
            ris(3) = ris(3) + x1^2;
            
        end  
      
end %switch ind_RHS
 
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

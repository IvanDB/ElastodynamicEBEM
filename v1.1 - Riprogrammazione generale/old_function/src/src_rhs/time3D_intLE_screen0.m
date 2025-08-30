function riga = time3D_intLE_screen0(pb_param,x1,x2,thk)

%Inizializzazione del BLOCCO 
riga = mat2cell(zeros(3,24),3,3*ones(1,8)); 

%VELOCITÀ delle ONDE S
velC_S = pb_param.velS;

%VELOCITÀ delle ONDE P
velC_P = pb_param.velP;

%Definiamo le COMPONENTI del VETTORE DISTANZA r tra un PUNTO di CAMPO
%e un PUNTO SORGENTE fissato (nodo di Gauss-Hammer)
r1 = @(y1,y2) y1-x1;
r2 = @(y1,y2) y2-x2;

%Definiamo il MODULO del VETTORE DISTANZA r tra un PUNTO di CAMPO
%e un PUNTO SORGENTE fissato (nodo di Gauss-Hammer)
r = @(y1,y2) sqrt(r1(y1,y2).^2+r2(y1,y2).^2);

%Primo triangolo
Integrazione={0,1,@(y1) 1-y1,1};

%Secondo triangolo
Integrazione=vertcat(Integrazione,{0,1,0,@(y1) 1-y1});

%Terzo triangolo
Integrazione=vertcat(Integrazione,{0,1,@(y1) -1+y1,0});

%Quarto triangolo
Integrazione=vertcat(Integrazione,{0,1,-1,@(y1)-1+y1});

%Quinto triangolo
Integrazione=vertcat(Integrazione,{-1,0,@(y1) 1+y1,1});

%Sesto triangolo
Integrazione=vertcat(Integrazione,{-1,0,0,@(y1) 1+y1});

%Settimo triangolo
Integrazione=vertcat(Integrazione,{-1,0,@(y1) -1-y1,0});

%Ottavo triangolo
Integrazione=vertcat(Integrazione,{-1,0,-1,@(y1) -1-y1});

%Definizione delle funzioni I_*
funI_Pt1 = @(y1,y2) (heaviside(thk-r(y1,y2)/velC_P).*(thk-r(y1,y2)/velC_P).*(thk+r(y1,y2)/velC_P))./(2*r(y1,y2).^3);
funI_St1 = @(y1,y2) (heaviside(thk-r(y1,y2)/velC_S).*(thk-r(y1,y2)/velC_S).*(thk+r(y1,y2)/velC_S))./(2*r(y1,y2).^3);

%Definizione delle funzioni J_*
funJ_St1 = @(y1,y2) heaviside(thk-r(y1,y2)/velC_S)./(velC_S^2*r(y1,y2));
funJ_Pt1 = @(y1,y2) heaviside(thk-r(y1,y2)/velC_P)./(velC_P^2*r(y1,y2));

%Definizione delle funzioni tilde_I_*
funtI_Pt1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funI_Pt1(y1,y2)./(r(y1,y2).^2);
funtI_St1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funI_St1(y1,y2)./(r(y1,y2).^2);
funtI_Pt1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_Pt1(y1,y2)./(r(y1,y2).^2);
funtI_St1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funI_St1(y1,y2)./(r(y1,y2).^2);
funtI_Pt1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funI_Pt1(y1,y2)./(r(y1,y2).^2);
funtI_St1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funI_St1(y1,y2)./(r(y1,y2).^2);


%Definizione delle funzioni tilde_J_*
funtJ_Pt1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funJ_Pt1(y1,y2)./(r(y1,y2).^2);
funtJ_St1_1_1 = @(y1,y2) (r1(y1,y2).^2).*funJ_St1(y1,y2)./(r(y1,y2).^2);
funtJ_Pt1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_Pt1(y1,y2)./(r(y1,y2).^2);
funtJ_St1_1_2 = @(y1,y2) (r1(y1,y2).*r2(y1,y2)).*funJ_St1(y1,y2)./(r(y1,y2).^2);
funtJ_Pt1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funJ_Pt1(y1,y2)./(r(y1,y2).^2);
funtJ_St1_2_2 = @(y1,y2) (r2(y1,y2).^2).*funJ_St1(y1,y2)./(r(y1,y2).^2);


nu_1_1=@(y1,y2)(funtJ_Pt1_1_1(y1,y2)-funtJ_St1_1_1(y1,y2))-(funI_Pt1(y1,y2)-funI_St1(y1,y2))...
            +3*(funtI_Pt1_1_1(y1,y2)-funtI_St1_1_1(y1,y2))+funJ_St1(y1,y2);
        
nu_2_2=@(y1,y2)(funtJ_Pt1_2_2(y1,y2)-funtJ_St1_2_2(y1,y2))-(funI_Pt1(y1,y2)-funI_St1(y1,y2))...
            +3*(funtI_Pt1_2_2(y1,y2)-funtI_St1_2_2(y1,y2))+funJ_St1(y1,y2);

nu_1_2=@(y1,y2)(funtJ_Pt1_1_2(y1,y2)-funtJ_St1_1_2(y1,y2))+3*(funtI_Pt1_1_2(y1,y2)-funtI_St1_1_2(y1,y2));

nu_3_3=@(y1,y2)-(funI_Pt1(y1,y2)-funI_St1(y1,y2))+funJ_St1(y1,y2);

   
for i=1:8
    ris=zeros(3,3);
    
    y1_min=Integrazione{i,1};
    y1_max=Integrazione{i,2};
    y2_min=Integrazione{i,3};
    y2_max=Integrazione{i,4};
    
    ris(1,1)=integral2(nu_1_1,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(1,1)=ris(1,1)/(4*pi*pb_param.rho);
    
    ris(2,2)=integral2(nu_2_2,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(2,2)=ris(2,2)/(4*pi*pb_param.rho);
    %N.B. Senza 'method','iterated' oppure 'AbsTol',1e-5,'RelTol',1e-6
    %compaiono dei NaN nel risultato dell'integrale ris(1,2)
    ris(1,2)=integral2(nu_1_2,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,2)=ris(1,2)/(4*pi*pb_param.rho);
    
    ris(2,1)=ris(1,2);
    
    ris(3,3)=integral2(nu_3_3,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(3,3)=ris(3,3)/(4*pi*pb_param.rho);
    
    riga{1,i}=ris;
end

riga=cell2mat(riga);

return
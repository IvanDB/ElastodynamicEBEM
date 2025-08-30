function riga = time3D_intLE_bar1(pb_param,x1,x2,x3,thk)

%Inizializzazione del BLOCCO 
riga = mat2cell(zeros(3,96),3,3*ones(1,32)); 

%VELOCITÀ delle ONDE S
velC_S = pb_param.velS;

%VELOCITÀ delle ONDE P
velC_P = pb_param.velP;

%Definiamo le COMPONENTI del VETTORE DISTANZA r tra un PUNTO di CAMPO
%e un PUNTO SORGENTE fissato (nodo di Gauss-Hammer)
r1 = @(y1,y2,y3) y1-x1;
r2 = @(y1,y2,y3) y2-x2;
r3 = @(y1,y2,y3) y3-x3;

%Definiamo il MODULO del VETTORE DISTANZA r tra un PUNTO di CAMPO
%e un PUNTO SORGENTE fissato (nodo di Gauss-Hammer)
r = @(y1,y2,y3) sqrt(r1(y1,y2,y3).^2+r2(y1,y2,y3).^2+r3(y1,y2,y3).^2);

%Definizione delle funzioni I_*
funI_Pt1 = @(y1,y2,y3) (heaviside(thk-r(y1,y2,y3)/velC_P).*(thk-r(y1,y2,y3)/velC_P).*(thk+r(y1,y2,y3)/velC_P))./(2*r(y1,y2,y3).^3);
funI_St1 = @(y1,y2,y3) (heaviside(thk-r(y1,y2,y3)/velC_S).*(thk-r(y1,y2,y3)/velC_S).*(thk+r(y1,y2,y3)/velC_S))./(2*r(y1,y2,y3).^3);

%Definizione delle funzioni J_*
funJ_St1 = @(y1,y2,y3) heaviside(thk-r(y1,y2,y3)/velC_S)./(velC_S^2*r(y1,y2,y3));
funJ_Pt1 = @(y1,y2,y3) heaviside(thk-r(y1,y2,y3)/velC_P)./(velC_P^2*r(y1,y2,y3));

%Definizione delle funzioni tilde_I_*
funtI_Pt1_1_1 = @(y1,y2,y3) (r1(y1,y2,y3).^2).*funI_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtI_St1_1_1 = @(y1,y2,y3) (r1(y1,y2,y3).^2).*funI_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

funtI_Pt1_2_2 = @(y1,y2,y3) (r2(y1,y2,y3).^2).*funI_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtI_St1_2_2 = @(y1,y2,y3) (r2(y1,y2,y3).^2).*funI_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

funtI_Pt1_3_3 = @(y1,y2,y3) (r3(y1,y2,y3).^2).*funI_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtI_St1_3_3 = @(y1,y2,y3) (r3(y1,y2,y3).^2).*funI_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

funtI_Pt1_1_2 = @(y1,y2,y3) (r1(y1,y2,y3).*r2(y1,y2,y3)).*funI_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtI_St1_1_2 = @(y1,y2,y3) (r1(y1,y2,y3).*r2(y1,y2,y3)).*funI_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

funtI_Pt1_1_3 = @(y1,y2,y3) (r1(y1,y2,y3).*r3(y1,y2,y3)).*funI_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtI_St1_1_3 = @(y1,y2,y3) (r1(y1,y2,y3).*r3(y1,y2,y3)).*funI_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

funtI_Pt1_2_3 = @(y1,y2,y3) (r2(y1,y2,y3).*r3(y1,y2,y3)).*funI_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtI_St1_2_3 = @(y1,y2,y3) (r2(y1,y2,y3).*r3(y1,y2,y3)).*funI_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

%Definizione delle funzioni tilde_J_*
funtJ_Pt1_1_1 = @(y1,y2,y3) (r1(y1,y2,y3).^2).*funJ_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtJ_St1_1_1 = @(y1,y2,y3) (r1(y1,y2,y3).^2).*funJ_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

funtJ_Pt1_2_2 = @(y1,y2,y3) (r2(y1,y2,y3).^2).*funJ_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtJ_St1_2_2 = @(y1,y2,y3) (r2(y1,y2,y3).^2).*funJ_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

funtJ_Pt1_3_3 = @(y1,y2,y3) (r3(y1,y2,y3).^2).*funJ_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtJ_St1_3_3 = @(y1,y2,y3) (r3(y1,y2,y3).^2).*funJ_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

funtJ_Pt1_1_2 = @(y1,y2,y3) (r1(y1,y2,y3).*r2(y1,y2,y3)).*funJ_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtJ_St1_1_2 = @(y1,y2,y3) (r1(y1,y2,y3).*r2(y1,y2,y3)).*funJ_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

funtJ_Pt1_1_3 = @(y1,y2,y3) (r1(y1,y2,y3).*r3(y1,y2,y3)).*funJ_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtJ_St1_1_3 = @(y1,y2,y3) (r1(y1,y2,y3).*r3(y1,y2,y3)).*funJ_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

funtJ_Pt1_2_3 = @(y1,y2,y3) (r2(y1,y2,y3).*r3(y1,y2,y3)).*funJ_Pt1(y1,y2,y3)./(r(y1,y2,y3).^2);
funtJ_St1_2_3 = @(y1,y2,y3) (r2(y1,y2,y3).*r3(y1,y2,y3)).*funJ_St1(y1,y2,y3)./(r(y1,y2,y3).^2);

%Definizione delle funzioni integrande
nu_1_1=@(y1,y2,y3)(funtJ_Pt1_1_1(y1,y2,y3)-funtJ_St1_1_1(y1,y2,y3))-(funI_Pt1(y1,y2,y3)-funI_St1(y1,y2,y3))...
            +3*(funtI_Pt1_1_1(y1,y2,y3)-funtI_St1_1_1(y1,y2,y3))+funJ_St1(y1,y2,y3);
        
nu_2_2=@(y1,y2,y3)(funtJ_Pt1_2_2(y1,y2,y3)-funtJ_St1_2_2(y1,y2,y3))-(funI_Pt1(y1,y2,y3)-funI_St1(y1,y2,y3))...
            +3*(funtI_Pt1_2_2(y1,y2,y3)-funtI_St1_2_2(y1,y2,y3))+funJ_St1(y1,y2,y3);

nu_3_3=@(y1,y2,y3)(funtJ_Pt1_3_3(y1,y2,y3)-funtJ_St1_3_3(y1,y2,y3))-(funI_Pt1(y1,y2,y3)-funI_St1(y1,y2,y3))...
            +3*(funtI_Pt1_3_3(y1,y2,y3)-funtI_St1_3_3(y1,y2,y3))+funJ_St1(y1,y2,y3);

nu_1_2=@(y1,y2,y3)(funtJ_Pt1_1_2(y1,y2,y3)-funtJ_St1_1_2(y1,y2,y3))+3*(funtI_Pt1_1_2(y1,y2,y3)-funtI_St1_1_2(y1,y2,y3));

nu_1_3=@(y1,y2,y3)(funtJ_Pt1_1_3(y1,y2,y3)-funtJ_St1_1_3(y1,y2,y3))+3*(funtI_Pt1_1_3(y1,y2,y3)-funtI_St1_1_3(y1,y2,y3));

nu_2_3=@(y1,y2,y3)(funtJ_Pt1_2_3(y1,y2,y3)-funtJ_St1_2_3(y1,y2,y3))+3*(funtI_Pt1_2_3(y1,y2,y3)-funtI_St1_2_3(y1,y2,y3));


%% Faccia superiore nel piano z=1

%-----------------------------------------------------------------------
%1° triangolo
Integrazione={0,1,-1,@(y1,y2,y3)-1+y1,1,1};
%2° triangolo
Integrazione=vertcat(Integrazione,{0,1,@(y1,y2)-1+y1,0,1,1});
%3° triangolo
Integrazione=vertcat(Integrazione,{0,1,0,@(y1,y2)y1,1,1});
%4° triangolo
Integrazione=vertcat(Integrazione,{0,1,@(y1,y2)y1,1,1,1});
%5° triangolo
Integrazione=vertcat(Integrazione,{-1,0,-1,@(y1,y2)y1,1,1});
%6° triangolo
Integrazione=vertcat(Integrazione,{-1,0,@(y1,y2)y1,0,1,1});
%7° triangolo
Integrazione=vertcat(Integrazione,{-1,0,0,@(y1,y2) 1+y1,1,1});
%8° triangolo
Integrazione=vertcat(Integrazione,{-1,0,@(y1,y2)1+y1,1,1,1});
%-----------------------------------------------------------------------
for i=1:8
       
    ris=zeros(3,3);
    
    y1_min=Integrazione{i,1};
    y1_max=Integrazione{i,2};
    y2_min=Integrazione{i,3};
    y2_max=Integrazione{i,4};
    y3_min=Integrazione{i,5};
    
    nu_11=@(y1,y2)nu_1_1(y1,y2,y3_min);
    nu_22=@(y1,y2)nu_2_2(y1,y2,y3_min);
    nu_33=@(y1,y2)nu_3_3(y1,y2,y3_min);
    nu_12=@(y1,y2)nu_1_2(y1,y2,y3_min);
    nu_13=@(y1,y2)nu_1_3(y1,y2,y3_min);
    nu_23=@(y1,y2)nu_2_3(y1,y2,y3_min);
    
    ris(1,1)=integral2(nu_11,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(1,1)=ris(1,1)/(4*pi*pb_param.rho);
    
    ris(2,2)=integral2(nu_22,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(2,2)=ris(2,2)/(4*pi*pb_param.rho);
    
    ris(3,3)=integral2(nu_33,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(3,3)=ris(3,3)/(4*pi*pb_param.rho);

    %N.B. Senza 'method','iterated' oppure 'AbsTol',1e-5,'RelTol',1e-6
    %compaiono dei NaN nel risultato dell'integrale ris(1,2)
    ris(1,2)=integral2(nu_12,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,2)=ris(1,2)/(4*pi*pb_param.rho);
    ris(2,1)=ris(1,2);
    
    ris(1,3)=integral2(nu_13,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,3)=ris(1,3)/(4*pi*pb_param.rho);
    ris(3,1)=ris(1,3);
    
    ris(2,3)=integral2(nu_23,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(2,3)=ris(2,3)/(4*pi*pb_param.rho);
    ris(3,2)=ris(2,3);
       
    riga{1,i}=ris;
end

%% Faccia laterale nel piano x=1
%-----------------------------------------------------------------------
%9° triangolo
Integrazione=vertcat(Integrazione,{1,1,-1,0,0,@(y2,y3)-y2});
%10° triangolo
Integrazione=vertcat(Integrazione,{1,1,-1,0,@(y2,y3)-y2,1});
%11° triangolo
Integrazione=vertcat(Integrazione,{1,1,0,1,0,@(y2,y3)1-y2});
%12° triangolo
Integrazione=vertcat(Integrazione,{1,1,0,1,@(y2,y3)1-y2,1});
%-----------------------------------------------------------------------
for i=9:12
       
    ris=zeros(3,3);
    
    y1_min=Integrazione{i,1};
    %y1_max=Integrazione{i,2};
    y2_min=Integrazione{i,3};
    y2_max=Integrazione{i,4};
    y3_min=Integrazione{i,5};
    y3_max=Integrazione{i,6};
    
    nu_11=@(y2,y3)nu_1_1(y1_min,y2,y3);
    nu_22=@(y2,y3)nu_2_2(y1_min,y2,y3);
    nu_33=@(y2,y3)nu_3_3(y1_min,y2,y3);
    nu_12=@(y2,y3)nu_1_2(y1_min,y2,y3);
    nu_13=@(y2,y3)nu_1_3(y1_min,y2,y3);
    nu_23=@(y2,y3)nu_2_3(y1_min,y2,y3);
    
    ris(1,1)=integral2(nu_11,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(1,1)=ris(1,1)/(4*pi*pb_param.rho);
    
    ris(2,2)=integral2(nu_22,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(2,2)=ris(2,2)/(4*pi*pb_param.rho);
    
    ris(3,3)=integral2(nu_33,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(3,3)=ris(3,3)/(4*pi*pb_param.rho);

    %N.B. Senza 'method','iterated' oppure 'AbsTol',1e-5,'RelTol',1e-6
    %compaiono dei NaN nel risultato dell'integrale ris(1,2)
    ris(1,2)=integral2(nu_12,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,2)=ris(1,2)/(4*pi*pb_param.rho);
    ris(2,1)=ris(1,2);
    
    ris(1,3)=integral2(nu_13,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,3)=ris(1,3)/(4*pi*pb_param.rho);
    ris(3,1)=ris(1,3);
    
    ris(2,3)=integral2(nu_23,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(2,3)=ris(2,3)/(4*pi*pb_param.rho);
    ris(3,2)=ris(2,3);
       
    riga{1,i}=ris;
end

%% Faccia inferiore nel piano z=0

%-----------------------------------------------------------------------
%13° triangolo
Integrazione=vertcat(Integrazione,{0,1,-1,@(y1,y2)-1+y1,0,0});
%14° triangolo
Integrazione=vertcat(Integrazione,{0,1,@(y1,y2)-1+y1,0,0,0});
%15° triangolo
Integrazione=vertcat(Integrazione,{-1,0,-1,@(y1,y2)y1,0,0});
%16° triangolo
Integrazione=vertcat(Integrazione,{-1,0,@(y1,y2)y1,0,0,0});
%17° triangolo
Integrazione=vertcat(Integrazione,{0,1,0,@(y1,y2)y1,0,0});
%18° triangolo
Integrazione=vertcat(Integrazione,{0,1,@(y1,y2)y1,1,0,0});
%19° triangolo
Integrazione=vertcat(Integrazione,{-1,0,0,@(y1,y2) 1+y1,0,0});
%20° triangolo
Integrazione=vertcat(Integrazione,{-1,0,@(y1,y2)1+y1,1,0,0});
%-----------------------------------------------------------------------

for i=13:20
       
    ris=zeros(3,3);
    
    y1_min=Integrazione{i,1};
    y1_max=Integrazione{i,2};
    y2_min=Integrazione{i,3};
    y2_max=Integrazione{i,4};
    y3_min=Integrazione{i,5};
    
    nu_11=@(y1,y2)nu_1_1(y1,y2,y3_min);
    nu_22=@(y1,y2)nu_2_2(y1,y2,y3_min);
    nu_33=@(y1,y2)nu_3_3(y1,y2,y3_min);
    nu_12=@(y1,y2)nu_1_2(y1,y2,y3_min);
    nu_13=@(y1,y2)nu_1_3(y1,y2,y3_min);
    nu_23=@(y1,y2)nu_2_3(y1,y2,y3_min);
    
    ris(1,1)=integral2(nu_11,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(1,1)=ris(1,1)/(4*pi*pb_param.rho);
    
    ris(2,2)=integral2(nu_22,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(2,2)=ris(2,2)/(4*pi*pb_param.rho);
    
    ris(3,3)=integral2(nu_33,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(3,3)=ris(3,3)/(4*pi*pb_param.rho);

    %N.B. Senza 'method','iterated' oppure 'AbsTol',1e-5,'RelTol',1e-6
    %compaiono dei NaN nel risultato dell'integrale ris(1,2)
    ris(1,2)=integral2(nu_12,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,2)=ris(1,2)/(4*pi*pb_param.rho);
    ris(2,1)=ris(1,2);
    
    ris(1,3)=integral2(nu_13,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,3)=ris(1,3)/(4*pi*pb_param.rho);
    ris(3,1)=ris(1,3);
    
    ris(2,3)=integral2(nu_23,y1_min,y1_max,y2_min,y2_max,'method','iterated','AbsTol',1e-5,'RelTol',1e-6);
    ris(2,3)=ris(2,3)/(4*pi*pb_param.rho);
    ris(3,2)=ris(2,3);
       
    riga{1,i}=ris;
end



%% Faccia laterale nel piano x=-1

%-----------------------------------------------------------------------
%21° triangolo
Integrazione=vertcat(Integrazione,{-1,-1,0,1,@(y2,y3)1-y2,1});
%22° triangolo
Integrazione=vertcat(Integrazione,{-1,-1,0,1,0,@(y2,y3)1-y2});
%23° triangolo
Integrazione=vertcat(Integrazione,{-1,-1,-1,0,@(y2,y3)-y2,1});
%24° triangolo
Integrazione=vertcat(Integrazione,{-1,-1,-1,0,0,@(y2,y3)-y2});
%-----------------------------------------------------------------------
for i=21:24
       
    ris=zeros(3,3);
    
    y1_min=Integrazione{i,1};
    %y1_max=Integrazione{i,2};
    y2_min=Integrazione{i,3};
    y2_max=Integrazione{i,4};
    y3_min=Integrazione{i,5};
    y3_max=Integrazione{i,6};
    
    nu_11=@(y2,y3)nu_1_1(y1_min,y2,y3);
    nu_22=@(y2,y3)nu_2_2(y1_min,y2,y3);
    nu_33=@(y2,y3)nu_3_3(y1_min,y2,y3);
    nu_12=@(y2,y3)nu_1_2(y1_min,y2,y3);
    nu_13=@(y2,y3)nu_1_3(y1_min,y2,y3);
    nu_23=@(y2,y3)nu_2_3(y1_min,y2,y3);
    
    ris(1,1)=integral2(nu_11,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(1,1)=ris(1,1)/(4*pi*pb_param.rho);
    
    ris(2,2)=integral2(nu_22,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(2,2)=ris(2,2)/(4*pi*pb_param.rho);
    
    ris(3,3)=integral2(nu_33,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(3,3)=ris(3,3)/(4*pi*pb_param.rho);

    %N.B. Senza 'method','iterated' oppure 'AbsTol',1e-5,'RelTol',1e-6
    %compaiono dei NaN nel risultato dell'integrale ris(1,2)
    ris(1,2)=integral2(nu_12,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,2)=ris(1,2)/(4*pi*pb_param.rho);
    ris(2,1)=ris(1,2);
    
    ris(1,3)=integral2(nu_13,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,3)=ris(1,3)/(4*pi*pb_param.rho);
    ris(3,1)=ris(1,3);
    
    ris(2,3)=integral2(nu_23,y2_min,y2_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(2,3)=ris(2,3)/(4*pi*pb_param.rho);
    ris(3,2)=ris(2,3);
       
    riga{1,i}=ris;
end


%% Faccia laterale nel piano y=1

%-----------------------------------------------------------------------
%25° triangolo
Integrazione=vertcat(Integrazione,{-1,0,1,1,@(y1,y3)1+y1,1});
%26° triangolo
Integrazione=vertcat(Integrazione,{-1,0,1,1,0,@(y1,y3)1+y1});
%27° triangolo
Integrazione=vertcat(Integrazione,{0,1,1,1,@(y1,y3)y1,1});
%28° triangolo
Integrazione=vertcat(Integrazione,{0,1,1,1,0,@(y1,y3)y1});
%-----------------------------------------------------------------------
for i=25:28
       
    ris=zeros(3,3);
    
    y1_min=Integrazione{i,1};
    y1_max=Integrazione{i,2};
    y2_min=Integrazione{i,3};
    %y2_max=Integrazione{i,4};
    y3_min=Integrazione{i,5};
    y3_max=Integrazione{i,6};
    
    nu_11=@(y1,y3)nu_1_1(y1,y2_min,y3);
    nu_22=@(y1,y3)nu_2_2(y1,y2_min,y3);
    nu_33=@(y1,y3)nu_3_3(y1,y2_min,y3);
    nu_12=@(y1,y3)nu_1_2(y1,y2_min,y3);
    nu_13=@(y1,y3)nu_1_3(y1,y2_min,y3);
    nu_23=@(y1,y3)nu_2_3(y1,y2_min,y3);
    
    ris(1,1)=integral2(nu_11,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(1,1)=ris(1,1)/(4*pi*pb_param.rho);
    
    ris(2,2)=integral2(nu_22,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(2,2)=ris(2,2)/(4*pi*pb_param.rho);
    
    ris(3,3)=integral2(nu_33,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(3,3)=ris(3,3)/(4*pi*pb_param.rho);

    %N.B. Senza 'method','iterated' oppure 'AbsTol',1e-5,'RelTol',1e-6
    %compaiono dei NaN nel risultato dell'integrale ris(1,2)
    ris(1,2)=integral2(nu_12,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,2)=ris(1,2)/(4*pi*pb_param.rho);
    ris(2,1)=ris(1,2);
    
    ris(1,3)=integral2(nu_13,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,3)=ris(1,3)/(4*pi*pb_param.rho);
    ris(3,1)=ris(1,3);
    
    ris(2,3)=integral2(nu_23,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(2,3)=ris(2,3)/(4*pi*pb_param.rho);
    ris(3,2)=ris(2,3);
       
    riga{1,i}=ris;
end



%% Faccia laterale nel piano y=-1

%-----------------------------------------------------------------------
%29° triangolo
Integrazione=vertcat(Integrazione,{-1,0,-1,-1,@(y1,y3)1+y1,1});
%30° triangolo
Integrazione=vertcat(Integrazione,{-1,0,-1,-1,0,@(y1,y3)1+y1});
%31° triangolo
Integrazione=vertcat(Integrazione,{0,1,-1,-1,@(y1,y3)y1,1});
%32° triangolo
Integrazione=vertcat(Integrazione,{0,1,-1,-1,0,@(y1,y3)y1});
%-----------------------------------------------------------------------
for i=29:32
       
    ris=zeros(3,3);
    
    y1_min=Integrazione{i,1};
    y1_max=Integrazione{i,2};
    y2_min=Integrazione{i,3};
    %y2_max=Integrazione{i,4};
    y3_min=Integrazione{i,5};
    y3_max=Integrazione{i,6};
    
    nu_11=@(y1,y3)nu_1_1(y1,y2_min,y3);
    nu_22=@(y1,y3)nu_2_2(y1,y2_min,y3);
    nu_33=@(y1,y3)nu_3_3(y1,y2_min,y3);
    nu_12=@(y1,y3)nu_1_2(y1,y2_min,y3);
    nu_13=@(y1,y3)nu_1_3(y1,y2_min,y3);
    nu_23=@(y1,y3)nu_2_3(y1,y2_min,y3);
    
    ris(1,1)=integral2(nu_11,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(1,1)=ris(1,1)/(4*pi*pb_param.rho);
    
    ris(2,2)=integral2(nu_22,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(2,2)=ris(2,2)/(4*pi*pb_param.rho);
    
    ris(3,3)=integral2(nu_33,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-5);
    ris(3,3)=ris(3,3)/(4*pi*pb_param.rho);

    %N.B. Senza 'method','iterated' oppure 'AbsTol',1e-5,'RelTol',1e-6
    %compaiono dei NaN nel risultato dell'integrale ris(1,2)
    ris(1,2)=integral2(nu_12,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,2)=ris(1,2)/(4*pi*pb_param.rho);
    ris(2,1)=ris(1,2);
    
    ris(1,3)=integral2(nu_13,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(1,3)=ris(1,3)/(4*pi*pb_param.rho);
    ris(3,1)=ris(1,3);
    
    ris(2,3)=integral2(nu_23,y1_min,y1_max,y3_min,y3_max,'method','iterated','AbsTol',1e-6,'RelTol',1e-6);
    ris(2,3)=ris(2,3)/(4*pi*pb_param.rho);
    ris(3,2)=ris(2,3);
       
    riga{1,i}=ris;
end

riga=cell2mat(riga);

return
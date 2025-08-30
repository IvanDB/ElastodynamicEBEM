function [rist5,rist6] = time3D_intT3EX_S_P(zeta,Theta1,Theta2,Info_tri2D,B)
%function [ris1t,ris2t,ris3t,ris4t] = time3D_intT3EX_P(zeta,p1,p2,B)

%La function time3D_intT3EX_S_P() prende in INPUT:
%- la variabile zeta contenente la distanza del punto sorgente dalla 
%  sua proiezione nel piano in cui giace il triangolo di campo
%
%- la variabile Theta1 che contiene l'ampiezza del primo angolo che 
%  individua il triangolo nel piano
%
%- la variabile Theta2 che contiene l'ampiezza del secondo angolo che 
%  individua il triangolo nel piano
%
%- la struct Info_Tri2D che contiene le informazioni geometriche sul
%  triangolo figlio corrente (lunghezza dei lati, ampiezza degli angoli,
%  valore del seno e del coseno degli angoli)
%
%- l'array 3D B di formato 3x3x6, costituito da 6 matrici 3x3 ciascuna 
%  delle quali contiene una delle 6 tipologie di coefficienti che
%  permettono di esprimere il prodotto r_ir_k al variare di i,k=1,2,3
%
%e restituisce in OUTPUT:
%- le matrici 3x3 rist5 e rist6 contenenti, rispettivamente, 
%  il valore degli integrali delle funzioni 1/r e r_ir_j/r^3 per i,j=1,2,3 
%  calcolati sul TRIANGOLO in presenza delle ONDE S e P
%
%-------------------------------------------------------------------------


if (zeta>1.0e-6)
    %In questo caso la distanza zeta tra il punto sorgente e il piano del
    %triangolo di campo è maggiore di zero e quindi il triangolo sorgente 
    %e il triangolo di campo NON sono COMPLANARI
    
    % STEP 1: CALCOLIAMO I COEFFICIENTI nel caso zeta>0 

    %Calcoliamo i COEFFICIENTI G_i utili a calcolare gli INTEGRALI 
    %ANALITICI sul TRIANGOLO in presenza delle ONDE S e P nel caso
    %zeta>0(vedi Milroy)
    G = time3D_coeffG_T3_zeta_maggiore_0(zeta,Theta1,Theta2,Info_tri2D);
    %---------------------------------------------------------------------

    % STEP 2: CALCOLIAMO GLI INTEGRALI sul TRIANGOLO in presenza delle 
    %ONDE S e P nel caso zeta>0

    %INTEGRALE di 1/r sul TRIANGOLO in presenza delle ONDE S e P nel caso
    %zeta>0
    rist5 = G(1)*eye(3,3);   
%     R=@(theta)(Info_tri2D.a*Info_tri2D.sin_beta)./(sin(theta+Info_tri2D.beta));
%     f5=@(theta,rho) rho./(rho.^2+zeta^2).^(1/2);
%     RIS5=integral2(f5,Theta1,Theta2,0,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*eye(3,3);
%     diff5=abs(RIS5-rist5)

    %INTEGRALE di r_i*r_j/r^3 sul TRIANGOLO in presenza delle ONDE S e P
    %nel caso zeta>0
    rist6 = G(4)*B(:,:,1)+G(5)*B(:,:,2)+G(6)*B(:,:,3)+G(7)*B(:,:,4)+...
        +G(8)*B(:,:,5)+G(9)*B(:,:,6);    
%     f6=@(theta,rho) rho./(rho.^2+zeta^2).^(3/2);
%     RIS6=integral2(f6,Theta1,Theta2,0,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,1);
%     f6=@(theta,rho) (sin(Info_tri2D.gamma-theta).*rho.^2)./(rho.^2+zeta^2).^(3/2);
%     RIS6=RIS6+integral2(f6,Theta1,Theta2,0,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,2);
%     f6=@(theta,rho) (sin(theta).*rho.^2)./(rho.^2+zeta^2).^(3/2);
%     RIS6=RIS6+integral2(f6,Theta1,Theta2,0,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,3);
%     f6=@(theta,rho) ((sin(Info_tri2D.gamma-theta).^2).*rho.^3)./(rho.^2+zeta^2).^(3/2);
%     RIS6=RIS6+integral2(f6,Theta1,Theta2,0,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,4);
%     f6=@(theta,rho) ((sin(Info_tri2D.gamma-theta).*sin(theta)).*rho.^3)./(rho.^2+zeta^2).^(3/2);
%     RIS6=RIS6+integral2(f6,Theta1,Theta2,0,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,5);
%     f6=@(theta,rho) ((sin(theta).^2).*rho.^3)./(rho.^2+zeta^2).^(3/2);
%     RIS6=RIS6+integral2(f6,Theta1,Theta2,0,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,6);
%     diff6=abs(RIS6-rist6)

    %---------------------------------------------------------------------
 
else
    %In questo caso la distanza zeta tra il punto sorgente e il piano del
    %triangolo di campo è uguale a zero e quindi il triangolo sorgente 
    %e il triangolo di campo sono COMPLANARI
    
    % STEP 1: CALCOLIAMO I COEFFICIENTI nel caso zeta=0

    %Calcoliamo i COEFFICIENTI G_i utili a calcolare gli INTEGRALI 
    %ANALITICI sul TRIANGOLO in presenza delle ONDE S e P nel caso
    %zeta=0(vedi Milroy)
    G = time3D_coeffG_T3_zeta_uguale_0(Theta1,Theta2,Info_tri2D);
    %---------------------------------------------------------------------
  
    % STEP 2: CALCOLIAMO GLI INTEGRALI sul TRIANGOLO in presenza delle 
    %ONDE S e P nel caso zeta=0

    %INTEGRALE di 1/r sul TRIANGOLO in presenza delle ONDE S e P nel caso
    %zeta>0
    rist5 = G(1)*eye(3,3); 
    
%     F=Info_tri2D.a*Info_tri2D.sin_beta;
%     R=@(theta) F./sin(theta+Info_tri2D.beta);
%     f5=@(theta,rho)(rho./((rho.^2+zeta^2).^(1/2)));
%     RIST5=integral2(f5,Theta1,Theta2,0,R,'Method','iterated','AbsTol',1e-15,'RelTol',1e-15)*eye(3,3);
%     diff5=abs(RIST5-rist5)

    %INTEGRALE di r_i*r_j/r^3 sul TRIANGOLO in presenza delle ONDE S e P
    %nel caso zeta=0
    rist6 = G(7)*B(:,:,4)+G(8)*B(:,:,5)+G(9)*B(:,:,6); 
    %---------------------------------------------------------------------

end %Fine if (zeta>1.0e-6)


end
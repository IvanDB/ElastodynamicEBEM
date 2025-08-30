function [rist1,rist2,rist3,rist4] = time3D_intTCEX_S_P(zeta,R_S,pR_S,Theta1,Theta2,Info_tri2D,B)
%function [ris1t,ris2t,ris3t,ris4t] = time3D_intT3EX_P(zeta,p1,p2,B)

%La function time3D_intTCEX_S_P() prende in INPUT:
%- la variabile zeta contenente la distanza del punto sorgente dalla 
%  sua proiezione nel piano in cui giace il triangolo di campo
%
%- la variabile R_S=c_S*Delta contenente il raggio della superifice 
%  sferica individuata dal fronte d'onda relativo alle onde S
%
%- la variabile pR_S contenente il raggio del cerchio che rappresenta il 
%  fronte dell'onda S nel piano del triangolo di campo
%  pR_S = sqrt([c_S*Delta]^2-zeta^2)
%
%- la variabile Theta1 che contiene l'ampiezza del primo angolo che 
%  individua il "trapezio/triangolo circolare"
%
%- la variabile Theta2 che contiene l'ampiezza del secondo angolo che 
%  individua il "trapezio/triangolo circolare"
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
%- le matrici 3x3 ris1, ris2, ris3 e ris4 contenenti, rispettivamente, il
%  valore degli integrali delle funzioni 1/r, r_ir_j/r^3, 1/r^3 e
%  r_ir_j/r^5 per i,j=1,2,3 calcolati sul trapezio/triangolo circolare
%  in presenza della sola onda P 
%
%-------------------------------------------------------------------------


if (zeta>1.0e-6)
    % STEP 1: CALCOLIAMO I COEFFICIENTI nel caso zeta>0

    %Calcoliamo i COEFFICIENTI G_i utili a calcolare gli INTEGRALI 
    %ANALITICI SUL TRAPEZIO/TRAINGOLO CIRCOLARE in presenza delle ONDE S 
    %e P nel caso zeta>0 (vedi Milroy)
    G = time3D_coeffG_TC_zeta_maggiore_0(zeta,R_S,pR_S,Theta1,Theta2,Info_tri2D);
    
    %--------------------------------------------------------------------
    
    % STEP 2: CALCOLIAMO GLI INTEGRALI nel caso zeta>0

    %INTEGRALE di 1/r sul TRAPEZIO/TRIANGOLO CIRCOLARE in presenza della sola 
    %ONDA P nel caso zeta>0
    rist1 = G(1)*eye(3,3); 
%     R=@(theta)(Info_tri2D.a*Info_tri2D.sin_beta)./(sin(theta+Info_tri2D.beta));
%     f1=@(theta,rho) rho./(rho.^2+zeta^2).^(1/2);
%     RIS1=integral2(f1,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*eye(3,3);
%     diff1=abs(RIS1-rist1)


    %INTEGRALE di r_i*r_j/r^3 sul TRAPEZIO/TRIANGOLO CIRCOLARE in presenza 
    %della sola ONDA P nel caso zeta>0
    rist2 = G(4)*B(:,:,1)+G(5)*B(:,:,2)+G(6)*B(:,:,3)+G(7)*B(:,:,4)+...
            +G(8)*B(:,:,5)+G(9)*B(:,:,6);   
%     f2=@(theta,rho) rho./(rho.^2+zeta^2).^(3/2);
%     RIS2=integral2(f2,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,1);
%     f2=@(theta,rho) (sin(Info_tri2D.gamma-theta).*rho.^2)./(rho.^2+zeta^2).^(3/2);
%     RIS2=RIS2+integral2(f2,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,2);
%     f2=@(theta,rho) (sin(theta).*rho.^2)./(rho.^2+zeta^2).^(3/2);
%     RIS2=RIS2+integral2(f2,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,3);
%     f2=@(theta,rho) ((sin(Info_tri2D.gamma-theta).^2).*rho.^3)./(rho.^2+zeta^2).^(3/2);
%     RIS2=RIS2+integral2(f2,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,4);
%     f2=@(theta,rho) ((sin(Info_tri2D.gamma-theta).*sin(theta)).*rho.^3)./(rho.^2+zeta^2).^(3/2);
%     RIS2=RIS2+integral2(f2,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,5);
%     f2=@(theta,rho) ((sin(theta).^2).*rho.^3)./(rho.^2+zeta^2).^(3/2);
%     RIS2=RIS2+integral2(f2,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,6);
%     diff2=abs(RIS2-rist2)

    %INTEGRALE di 1/r^3 sul sul TRAPEZIO/TRIANGOLO CIRCOLARE in presenza 
    %della sola ONDA P nel caso zeta>0
    rist3 = G(4)*eye(3,3); 
%     R=@(theta)(Info_tri2D.a*Info_tri2D.sin_beta)./(sin(theta+Info_tri2D.beta));
%     f3=@(theta,rho) rho./(rho.^2+zeta^2).^(3/2);
%     RIS3=integral2(f3,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*eye(3,3);
%     diff3=abs(RIS3-rist3)

    %INTEGRALE di r_i*r_j/r^5 sul TRAPEZIO/TRIANGOLO CIRCOLARE in presenza 
    %della sola ONDA P nel caso zeta>0
    rist4 = G(14)*B(:,:,1)+G(15)*B(:,:,2)+G(16)*B(:,:,3)+G(17)*B(:,:,4)+...
        +G(18)*B(:,:,5)+G(19)*B(:,:,6);   
%     f4=@(theta,rho) rho./(rho.^2+zeta^2).^(5/2);
%     RIS4=integral2(f4,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,1);
%     f4=@(theta,rho) (sin(Info_tri2D.gamma-theta).*rho.^2)./(rho.^2+zeta^2).^(5/2);
%     RIS4=RIS4+integral2(f4,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,2);
%     f4=@(theta,rho) (sin(theta).*rho.^2)./(rho.^2+zeta^2).^(5/2);
%     RIS4=RIS4+integral2(f4,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,3);
%     f4=@(theta,rho) ((sin(Info_tri2D.gamma-theta).^2).*rho.^3)./(rho.^2+zeta^2).^(5/2);
%     RIS4=RIS4+integral2(f4,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,4);
%     f4=@(theta,rho) ((sin(Info_tri2D.gamma-theta).*sin(theta)).*rho.^3)./(rho.^2+zeta^2).^(5/2);
%     RIS4=RIS4+integral2(f4,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,5);
%     f4=@(theta,rho) ((sin(theta).^2).*rho.^3)./(rho.^2+zeta^2).^(5/2);
%     RIS4=RIS4+integral2(f4,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*B(:,:,6);
%     diff4=abs(RIS4-rist4)

    
    %--------------------------------------------------------------------
  
else
    % STEP 1: CALCOLIAMO I COEFFICIENTI nel caso zeta=0

    %Calcoliamo i COEFFICIENTI G_i utili a calcolare gli INTEGRALI 
    %ANALITICI SUL TRAPEZIO/TRAINGOLO CIRCOLARE in presenza delle ONDE S 
    %e P nel caso zeta=0 (vedi Milroy)
    G = time3D_coeffG_TC_zeta_uguale_0(R_S,Theta1,Theta2,Info_tri2D);
    
    %--------------------------------------------------------------------
    
    % STEP 2: CALCOLIAMO GLI INTEGRALI nel caso zeta=0

    %INTEGRALE di 1/r sul TRAPEZIO/TRIANGOLO CIRCOLARE in presenza della sola 
    %ONDA P nel caso zeta=0
    rist1 = G(1)*eye(3,3); 

%     R=@(theta)(Info_tri2D.a*Info_tri2D.sin_beta)./(sin(theta+Info_tri2D.beta));
%     f1=@(theta,rho) rho./(rho.^2+zeta^2).^(1/2);
%     RIS1=integral2(f1,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*eye(3,3);
%     diff1=abs(RIS1-rist1)


    %INTEGRALE di r_i*r_j/r^3 sul TRAPEZIO/TRIANGOLO CIRCOLARE in presenza 
    %della sola ONDA P nel caso zeta=0
    rist2 = G(7)*B(:,:,4)+G(8)*B(:,:,5)+G(9)*B(:,:,6);

    %INTEGRALE di 1/r^3 sul sul TRAPEZIO/TRIANGOLO CIRCOLARE in presenza 
    %della sola ONDA P nel caso zeta=0
    rist3 = G(4)*eye(3,3); 

%     R=@(theta)(Info_tri2D.a*Info_tri2D.sin_beta)./(sin(theta+Info_tri2D.beta));
%     f3=@(theta,rho) rho./(rho.^2+zeta^2).^(3/2);
%     RIS3=integral2(f3,Theta1,Theta2,pR_S,R,'Method','iterated','AbsTol',0,'RelTol',1e-15)*eye(3,3);
%     diff3=abs(RIS3-rist3)

    %INTEGRALE di r_i*r_j/r^5 sul TRAPEZIO/TRIANGOLO CIRCOLARE in presenza 
    %della sola ONDA P nel caso zeta=0
    rist4 = G(17)*B(:,:,4)+G(18)*B(:,:,5)+G(19)*B(:,:,6);
    
    %--------------------------------------------------------------------

end %Fine if(zeta>1.0e-6)

return
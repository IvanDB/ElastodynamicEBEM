function [RIST]=rotazione_inversa_attorno_asse_z(rist,theta0)

%La function rotazione_inversa_attorno_asse_z() prende in input:
%- la matrice 3x3 rist che contiene i risultati degli integrali di
%  r_ir_j/r^3 e r_ir_j/r^5 calcolati sul triangolo figlio ruotato 
%  che vanno moltiplicati per i coefficienti di "correzione" derivanti 
%  dall'applicazione della rotazione, in modo tale da ricondurre questi 
%  integrali ad essere integrali sul triangolo figlio non ruotato
%
%- la variabile theta0 contenente l'angolo (orientato in senso antiorario)
%  della rotazione diretta rispetto al terzo asse della terna cartesiana 
%  che porta il triangolo figlio nel triangolo figlio ruotato
%
%e restituisce in OUTPUT:
%- la matrice 3x3 RIST contenente i risultati degli integrali "corretti"
%  dalla rotazione attorno all'asse z nel sistema di riferimento del
%  triangolo di campo

%Calcoliamo il coseno e il seno dell'angolo di rotazione theta0
cos_theta0=cos(theta0);
sin_theta0=sin(theta0);

%Costruiamo la matrice A di rotazione inversa (rotazione in senso orario
%attorno all'asse z)
A=[cos_theta0 sin_theta0 0; -sin_theta0 cos_theta0 0; 0 0 1];

%Calcoliamo gli integrali sul triangolo figlio non ruotato attraverso
%questo prodotto di matrici
RIST=A*(rist*A');


%     RIST=zeros(3,3);
%     cos_theta0=cos(theta0);
%     sin_theta0=sin(theta0);
%     cos_theta0_2=cos_theta0^2;
%     sin_theta0_2=sin_theta0^2;
%     cos_sin_theta0=cos_theta0*sin_theta0;  
%     RIST(1,1)=rist(1,1)*cos_theta0_2+rist(2,2)*sin_theta0_2+2*rist(1,2)*cos_sin_theta0;
%     RIST(2,2)=rist(1,1)*sin_theta0_2+rist(2,2)*cos_theta0_2-2*rist(1,2)*cos_sin_theta0;
%     RIST(3,3)=rist(3,3); 
%     RIST(1,2)=(-rist(1,1)+rist(2,2))*cos_sin_theta0+rist(1,2)*(cos_theta0_2-sin_theta0_2);
%     RIST(2,1)=RIST(1,2);
%     RIST(1,3)=rist(1,3)*cos_theta0+rist(2,3)*sin_theta0;
%     RIST(3,1)=RIST(1,3); 
%     RIST(2,3)=-rist(1,3)*sin_theta0+rist(2,3)*cos_theta0;
%     RIST(3,2)= RIST(2,3);
    
%     RIST(1,1)=rist(1,1)*cos_theta0_2+rist(2,2)*sin_theta0_2+rist(1,2)*cos_sin_theta0+rist(2,1)*cos_sin_theta0;
%     RIST(2,2)=rist(1,1)*sin_theta0_2+rist(2,2)*cos_theta0_2-rist(1,2)*cos_sin_theta0-rist(2,1)*cos_sin_theta0;
%     RIST(3,3)=rist(3,3);
%     
%     RIST(1,2)=(-rist(1,1)+rist(2,2))*cos_sin_theta0+rist(1,2)*cos_theta0_2-rist(2,1)*sin_theta0_2;
%     RIST(2,1)=(-rist(1,1)+rist(2,2))*cos_sin_theta0+rist(2,1)*cos_theta0_2-rist(1,2)*sin_theta0_2;
%     
%     RIST(1,3)=rist(1,3)*cos_theta0+rist(2,3)*sin_theta0;
%     RIST(3,1)=rist(3,1)*cos_theta0+rist(3,2)*sin_theta0;
%     
%     RIST(2,3)=-rist(1,3)*sin_theta0+rist(2,3)*cos_theta0;
%     RIST(3,2)=-rist(3,1)*sin_theta0+rist(3,2)*cos_theta0;
return
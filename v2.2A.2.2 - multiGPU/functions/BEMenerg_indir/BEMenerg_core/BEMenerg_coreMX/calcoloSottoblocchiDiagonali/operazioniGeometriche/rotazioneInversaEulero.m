function ris = rotazioneInversaEulero(SdR, vnF, ris)
% INPUT
% - 
% - 
% OUTPUT
% - 

%% PREPARAZIONE VERSORI

%Estrazione versori terna relativa al triangolo di campo
X = SdR.e1;
Z = vnF;

%Defizione versori SdR canonico di R3
e1 = [1 0 0];
e3 = [0 0 1];

%% CALCOLO ANGOLI DI EULERO

%Calcolo LINEA dei NODI

%Check buona definizione come intersezione di piani
if (norm(cross(e3, Z), 2) < 1.0e-10)
    %Caso di piani XY ed e1e2 paralleli

    %Convenzione sulla scelta della linea
    l_nodi = X;    
    
    %Calcolo angolo di precessione phi
    cos_phi = e1 * l_nodi';
    phi = acos(cos_phi);
    
    %Check per eventuale aggiustamento di phi dovuto al codominio di acos
    if (cross(e1, l_nodi) * e3' < -1.0e-10)
        phi = 2*pi - phi;
    end

    %Calcolo angolo di nutazione (caso banale (piani paralleli))
    omega = pi * (e3*Z'< 0);

    %Calcolo angolo di rotazione propria (caso banale (piani paralleli))
    psi = 0;

else
    %Caso di piani XY ed e1e2 non paralleli

    l_nodi = cross(e3, Z) / norm(cross(e3, Z), 2);   
    
    %Calcolo angolo di precessione phi
    cos_phi = e1 * l_nodi';
    phi = acos(cos_phi);
    
    %Check per eventuale aggiustamento di phi dovuto al codominio di acos
    if (cross(e1, l_nodi) * e3' < -1.0e-10)
        phi = 2*pi - phi;
    end

    %Calcolo angolo di nutazione (caso non banale)
    cos_omega = e3 * Z';
    omega = acos(cos_omega);

    
    %Calcolo angolo di rotazione propria (caso non banale)
    cos_psi = l_nodi * X';
    psi = acos(cos_psi);
    psi = real(psi);

    %Check per eventuale aggiustamento di phi dovuto al codominio di acos
    if (cross(l_nodi, X) * Z'< -1.0e-10)
        psi = 2*pi - psi;
    end

end

%% COSTRUZIONE MATRICI DI ROTAZIONE

A_psi   = [cos(psi) sin(psi)    0;
          -sin(psi) cos(psi)    0;
           0        0           1];
A_omega = [1        0           0;
           0        cos(omega)  sin(omega); 
           0        -sin(omega) cos(omega)];
A_phi   = [cos(phi) sin(phi)    0;
          -sin(phi) cos(phi)    0; 
           0        0           1];

%Calcolo matrice di rotazione complessiva
C = A_psi * A_omega * A_phi;

%Pulizia numerica della matrice di rotazione
C(abs(C) < 1.0e-14) = 0;

%% APPLICAZIONE MATRICE DI ROTAZIONE
ris = C' * (ris * C);
return
end
function dist = calcoloDistanzaTriangolo(tri2D)

%Calcolo DISTANZA VERTICI TRIANGOLO di CAMPO dall'ORIGINE (0,0)
D1 = sqrt(sum(tri2D(1, :).^2));
D2 = sqrt(sum(tri2D(2, :).^2));
D3 = sqrt(sum(tri2D(3, :).^2));  

%Inizializzazione matrice 3x6 contenente in ogni riga distanze e coordinate
% di ciascuna delle 3 coppie di vertici possibili 
A = [D1, D2, tri2D(1, :), tri2D(2, :); 
     D2, D3, tri2D(2, :), tri2D(3, :);
     D1, D3, tri2D(1, :), tri2D(3, :)];

%Inizializzazione vettore contente le distanze di ciascun lato
d_min = zeros(3, 1);

%Ciclo sui 3 lati del triangolo
for i = 1 : 3
    
    %Estrazione coordinate coppia vertici del lato corrente
    P1 = A(i, 3:4);
    P2 = A(i, 5:6);
    
    %Estrazione distanze coppia vertici del lato corrente
    d1 = A(i, 1);
    d2 = A(i, 2);

    %Calcolo vettore congiungente la coppia di vertici corrente
    vL = P2 - P1;
    %Calcolo lunghezza vettore congiungente la coppia di vertici corrente
    L = sqrt(sum(vL.^2));

    %Calcolo distanza dell'origine dal segmento p1p2
    if (-P1*vL' < -1.0e-6)
        %Caso angolo OP1P2 ottuso
        d_min(i) = d1;

    elseif (-P2*(-vL)' < -1.0e-6)
        %Caso angolo OP2P1 ottuso
        d_min(i) = d2;

    else
        %Caso angoli OP1P2 e OP2P1 entrambi acuti

        cos_angl_OP1P2 = (-P1/d1) * (vL/L)';
        sin_angl_OP1P2 = sqrt(1 - cos_angl_OP1P2^2);
        d_min(i) = d1 * sin_angl_OP1P2;
    end 
    
    %Calcolo minima distanza complessiva
    dist = min(d_min);
end
return
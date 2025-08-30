function dist=distanza_triangolo(tri2D)

%La function distanza_triangolo() prende in INPUT:
%- la matrice 3x2 tri2D in cui l'i-esima riga contiene le coordinate del 
%  vertice i-esimo del triangolo di campo nel nuovo sistema di 
%  riferimento
%
%e restituisce in OUTPUT:
%- la variabile dist contenente la minima distanza dell'origine (0,0) dal
%  triangolo di campo
%
%-------------------------------------------------------------------------
%N.B. Stiamo ragionando nell'ipotesi in cui l'orgine (0,0) NON appartiene
%al triangolo di campo
%-------------------------------------------------------------------------

%Calcoliamo la DISTANZA di ciascuno dei TRE VERTICI del 
%TRIANGOLO di CAMPO dall'ORIGINE (0,0) del nuovo SISTEMA di RIFERIMENTO 
%BIDIMENSIONALE 
D1=sqrt(sum(tri2D(1,:).^2));
D2=sqrt(sum(tri2D(2,:).^2));
D3=sqrt(sum(tri2D(3,:).^2));  

%Creiamo una matrice 3x6 A in cui ogni riga fa riferimento ad una diversa
%coppia di vertici del triangolo di campo. I primi due elementi di 
%ciascuna riga corrispondono alla distanza dei due vertici considerati 
%dall'origine, mentre i restanti elementi contengono le coordinate dei 
%vertici considerati
A=[D1 D2 tri2D(1,:) tri2D(2,:); D2 D3 tri2D(2,:) tri2D(3,:); D1 D3 tri2D(1,:) tri2D(3,:)];

%Creiamo il vettore 3x1 d_min in cui memorizzeremo la minima distanza
%dell'origine (0,0) da ciascuno dei tre lati del triangolo di campo
d_min=zeros(3,1);

for i=1:3
    %Nell'i-esima iterazione calcoliamo la distanza dell'orgine
    %dall'i-esimo lato del triangolo di campo
    
    %Memorizziamo nei vettori 1x2 p1 e p2 le coordinate dei due vertici 
    %del triangolo di campo che stiamo considerando
    p1=A(i,3:4);
    p2=A(i,5:6);
    
    %Memorizziamo nelle variabili d1 e d2 le distanze dall'origine dei 
    %due vertici del triangolo di campo che stiamo considerando
    d1=A(i,1);
    d2=A(i,2);

    %VETTORE DISTANZA p2-p1 (che è un vettore lungo uno dei lati del triangolo
    %di campo)
    vL = p2-p1;
    %LUNGHEZZA del VETTORE p2-p1
    L = sqrt(sum(vL.^2));

    %Calcoliamo la minima distanza dell'origine dal segmento che unisce i 
    %vertici p1 e p2 del triangolo
    if (-p1*vL'<-1.0e-6)

        %Se il prodotto scalalre tra -p1 e vL è negativo, significa che
        %l'angolo tra questi due vettori è ottuso. Quindi scegliamo come
        %distanza la distanza del vertice p1 dall'origine
        d_min(i)=d1;

    elseif (-p2*(-vL)'<-1.0e-6)
        %Se il prodotto scalalre tra -p2 e -vL è negativo, significa che
        %l'angolo tra questi due vettori è ottuso. Quindi scegliamo come
        %distanza la distanza del vertice p2 dall'origine
        d_min(i)=d2;

    else
        %Se ricadiamo in questo caso significa che nessuno degli angoli tra i
        %vettori -p1 e vL   e tra i vettori -p2 e -vL è ottuso. 
        %Quindi scegliamo la distanza dell'origine dalla retta passante per p1
        %e p2
        cos_angl=(-p1/d1)*(vL/L)';
        sin_angl=sqrt(1-cos_angl^2);
        d_min(i)=d1*sin_angl;
    end 
    
    %La  distanza minima dell'origine (0,0) dal triangolo di campo sarà
    %uguale al minimo tra le tre distanze dell'origine dai tre lati del
    %triangolo
    dist=min(d_min);
end


% %Creiamo una matrice 3x3 A in cui sulla riga i-esima sono presenti la
% %distanza dell'i-esimo vertice del triangolo di campo dall'origine e le  
% %sue due coordinate
% A=[d1 tri2D(1,:); d2 tri2D(2,:); d3 tri2D(3,:)];
% 
% %Ordiniamo le righe della matrice A rispetto alla prima colonna, ovvero
% %rispetto alle distanze dei tre vertici del triangolo di campo dall'origine
% A=sortrows(A);
% 
% %memorizziamo nei vettori 1x2 p1 e p2 le coordinate die due vertici del
% %triangolo di campo che hanno minore distanza dall'ogiine
% p1=A(1,2:end);
% p2=A(2,2:end);
% 
% 
% %VETTORE DISTANZA p2-p1 (che è un vettore lungo uno dei lati del triangolo
% %di campo)
% vL = p2-p1;
% %LUNGHEZZA del VETTORE p2-p1
% L = sqrt(sum(vL.^2));
% 
% %Calcoliamo la minima distanza dell'origine dal segmento che unisce i 
% %vertici p1 e p2 del triangolo
% if (-p1*vL'<-1.0e-6)
%     
%     %Se il prodotto scalalre tra -p1 e vL è negativo, significa che
%     %l'angolo tra questi due vettori è ottuso. Quindi scegliamo come
%     %distanza la distanza del vertice p1 dall'origine
%     dist=A(1,1);
%     
% elseif (-p2*(-vL)'<-1.0e-6)
%     %Se il prodotto scalalre tra -p2 e -vL è negativo, significa che
%     %l'angolo tra questi due vettori è ottuso. Quindi scegliamo come
%     %distanza la distanza del vertice p2 dall'origine
%     dist=A(2,2);
%     
% else
%     %Se ricadiamo in questo caso significa che nessuno degli angoli tra i
%     %vettori -p1 e vL   e tra i vettori -p2 e -vL è ottuso. 
%     %Quindi scegliamo la distanza dell'origine dalla retta passante per p1
%     %e p2
%     cos_angl=(-p1/A(1,1))*(vL/L)';
%     sin_angl=sqrt(1-cos_angl^2);
%     dist=A(1,1)*sin_angl;
% end  
   
return
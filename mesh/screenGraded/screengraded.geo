//******************* Definizione dei Parametri

lc =1;

//Lato del quadrato di base
L = 0.5;

//Numero di punti sul ciascun lato del quadrato [L,0]^2
N = 8;

B = 2;

//******************* Definizione dell'ostacolo per punti

//Primo lato

For k In {0:N}
	Point (k+1) = {L*(-1+(k/N)^B), -L, 0, lc};
EndFor

For k In {0:N-1}
	Point (2*N+1-k) = {L*(1-(k/N)^B), -L, 0, lc};
EndFor

//Secondo lato

For k In {1:N}
	Point (2*N+1+k) = {L, L*(-1+(k/N)^B), 0, lc};
EndFor

For k In {0:N-1}
	Point (4*N+1-k) = {L, L*(1-(k/N)^B), 0, lc};
EndFor

//Terzo lato

For k In {1:N}
	Point (4*N+1+k) = {L*(1-(k/N)^B), L, 0, lc};
EndFor

For k In {0:N-1}
	Point (6*N+1-k) = {L*(-1+(k/N)^B), L, 0, lc};
EndFor

//Quarto lato

For k In {1:N}
	Point (6*N+1+k) = {-L, L*(1-(k/N)^B), 0, lc};
EndFor

For k In {0:N-1}
	Point (8*N+1-k) = {-L, L*(-1+(k/N)^B), 0, lc};
EndFor

//******************* Definizione dell'ostacolo 

For k In {1:8*N-1}
	Line (k) = {k, k+1};
	coilLoops[k] = k;
EndFor

Line(8*N) = {8*N, 1};
coilLoops[8*N] = 8*N;

Line Loop(1) = {coilLoops[1]:coilLoops[8*N]};

Plane Surface(1)  = {1};
Surface Loop(2) = {1};

Transfinite Surface {1} = {1,2*N+1,4*N+1,6*N+1};

Physical Surface(1) = {1};

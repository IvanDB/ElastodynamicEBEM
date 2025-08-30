Mesh.ElementOrder = 1;

lc = 1.0;

//************************ Parametri del cubo

L = 1.0; //lunghezza
H = 1.0; //altezza
D = 0.0; //profondità

//************************ Punti base

Point(1) = {  L,  H, D,lc}; 
Point(2) = {  0,  H, D,lc}; 
Point(3) = { -L,  H, D,lc};  
Point(4) = { -L,  0, D,lc};  
Point(5) = { -L, -H, D,lc}; 
Point(6) = {  0, -H, D,lc};
Point(7) = {  L, -H, D,lc};  
Point(8) = {  L,  0, D,lc}; 


//************************ Linee base

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};

//************************ Superfici base

Line Loop(13) = {1,2,3,4,5,6,7,8};
Plane Surface(14) = {13};


Surface Loop(25) = {14};

//************************ Parametri di discretizzazione mesh base

//Transfinite Surface {14} = {7,5,3,1};

Physical Surface(1) = {14};
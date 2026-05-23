Mesh.ElementOrder = 1;

lc = 1.0;

//************************ Parametri del cubo

L = 0.5; //lunghezza
H = 0.5; //altezza
D = 1.0; //profondità

//************************ Punti base

Point(newp) = { L, H, D,lc};  /* Punto      1 */
Point(newp) = { L, H, 0,lc};  /* Punto      2 */
Point(newp) = {-L, H, D,lc};  /* Punto      3 */
Point(newp) = {-L,-H, D,lc};  /* Punto      4 */
Point(newp) = { L,-H, D,lc};  /* Punto      5 */
Point(newp) = { L,-H, 0,lc};  /* Punto      6 */
Point(newp) = {-L, H, 0,lc};  /* Punto      7 */
Point(newp) = {-L,-H, 0,lc};  /* Punto      8 */

Point(newp) = {0,0,0,lc};  /* Punto      9 */
Point(newp) = {0,0,D,lc};  /* Punto      10 */
Point(newp) = {0,H,D/2,lc};  /* Punto      11 */
Point(newp) = {0,-H,D/2,lc};  /* Punto      12 */

//************************ Linee base

Line(1) = {3,1};
Line(2) = {3,7};
Line(3) = {7,2};
Line(4) = {2,1};
Line(5) = {1,5};
Line(6) = {5,4};
Line(7) = {4,8};
Line(8) = {8,6};
Line(9) = {6,5};
Line(10) = {6,2};
Line(11) = {3,4};
Line(12) = {8,7};

//************************ Superfici base

Line Loop(13) = {-6,-5,-1,11};
Plane Surface(14) = {13};
Line Loop(15) = {4,5,-9,10};
Plane Surface(16) = {15};
Line Loop(17) = {3,12,-8,-10};
Plane Surface(18) = {17};
Line Loop(19) = {-7,-12,2,-11};
Plane Surface(20) = {19};
Line Loop(21) = {-4,-3,-2,1};
Plane Surface(22) = {21};
Line Loop(23) = {8,9,6,7};
Plane Surface(24) = {23};

Surface Loop(25) = {14,24,18,22,16,20};

//************************ Parametri di discretizzazione mesh base

Transfinite Surface {14} = {5,4,3,1};
Transfinite Surface {16} = {6,2,1,5};
Transfinite Surface {18} = {6,2,7,8};
Transfinite Surface {20} = {3,7,8,4};
Transfinite Surface {22} = {3,7,2,1};
Transfinite Surface {24} = {4,8,6,5};

Physical Surface(1) = {14,16,18,20,22,24};//+
Translate {0, 0, 1} {
  Point{3}; Curve{11}; Curve{6}; Surface{14}; Surface{24}; Point{4}; Point{8}; Point{12}; Point{10}; Point{5}; Point{7}; Point{9}; Point{6}; Point{11}; Point{1}; Point{2}; Curve{7}; Curve{12}; Curve{2}; Curve{9}; Curve{8}; Curve{5}; Curve{1}; Curve{4}; Curve{10}; Curve{3}; Surface{20}; Surface{22}; Surface{18}; Surface{16}; 
}

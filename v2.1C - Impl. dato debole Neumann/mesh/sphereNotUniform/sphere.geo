lc = 1.00;
Point (1) = { 0,   0 ,  0 , lc};Point (2) = {-1,   0 ,  0 , lc};Point (3) = { 1 ,  0 ,  0 , lc};Point (4) = { 0 , -1 ,  0 , lc};Point (5) = { 0 ,  1 ,  0 , lc};Point (6) = { 0 ,  0 ,  1 , lc};Point (7) = { 0 ,  0 , -1 , lc};

Circle ( 1) = {3 , 1 , 5};Circle ( 2) = {2 , 1 , 5};Circle ( 3) = {2 , 1 , 4};Circle ( 4) = {4 , 1 , 3};Circle ( 5) = {6 , 1 , 3};Circle ( 6) = {6 , 1 , 2};Circle ( 7) = {2 , 1 , 7};Circle ( 8) = {7 , 1 , 3};
Circle ( 9) = {5 , 1 , 6};
Circle (10) = {6 , 1 , 4};
Circle (11) = {4 , 1 , 7};
Circle (12) = {7 , 1 , 5};

Line Loop (13) = {9 , 5 , 1};
Ruled Surface (14) = {13};
Line Loop (15) = {2 , 9 , 6};
Ruled Surface (16) = {-15};Line Loop (17) = {7 , 12 , -2};
Ruled Surface (18) = {-17};
Line Loop (19) = {12 , -1 , -8};
Ruled Surface (20) = {19};
Line Loop (21) = {6 , 3 , -10};
Ruled Surface (22) = {21};
Line Loop (23) = {10 , 4 , -5};
Ruled Surface (24) = {23};
Line Loop (25) = {11 , -7 , 3};
Ruled Surface (26) = {-25};
Line Loop (27) = {8 , -4 , 11};
Ruled Surface (28) = {27};
Surface Loop (29) = {14 , 16 , 18 , 20 , 22 , 24 , 26 , 28};

Physical Surface(7) = {14 , 16 , 18 , 24 , 22 , 20 , 26 , 28};

//Volume (1) = {1};
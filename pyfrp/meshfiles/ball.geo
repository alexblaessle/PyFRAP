//Parameters
volSize_px=20.;                                                                   
radius=968.181242777;
center_x=256;
center_y=256;

//Points
Point(1) = {center_x, center_y, -radius, volSize_px};
Point(2) = {center_x+radius, center_y, -radius, volSize_px};
Point(3) = {center_x, center_y+radius, -radius, volSize_px};
Point(4) = {center_x-radius, center_y, -radius, volSize_px};
Point(5) = {center_x, center_y-radius, -radius, volSize_px};
Point(6) = {center_x, center_y, radius-radius, volSize_px};
Point(7) = {center_x, center_y, -radius-radius, volSize_px};

//Circles connecting points
Circle(1) = {3, 1, 7};
Circle(2) = {7, 1, 5};
Circle(3) = {5, 1, 6};
Circle(4) = {6, 1, 3};
Circle(5) = {3, 1, 4};
Circle(6) = {4, 1, 5};
Circle(7) = {5, 1, 2};
Circle(8) = {2, 1, 3};
Circle(9) = {2, 1, 7};
Circle(10) = {2, 1, 6};
Circle(11) = {6, 1, 4};
Circle(12) = {4, 1, 7};

//LineLoops and Surfaces
Line Loop(14) = {4, -8, 10};
Ruled Surface(14) = {14};
Line Loop(16) = {8, 1, -9};
Ruled Surface(16) = {16};
Line Loop(18) = {1, -12, -5};
Ruled Surface(18) = {18};
Line Loop(20) = {11, -5, -4};
Ruled Surface(20) = {20};
Line Loop(22) = {11, 6, 3};
Ruled Surface(22) = {22};
Line Loop(24) = {3, -10, -7};
Ruled Surface(24) = {24};
Line Loop(26) = {9, 2, 7};
Ruled Surface(26) = {26};
Line Loop(28) = {6, -2, -12};
Ruled Surface(28) = {28};

//Volume
Surface Loop(30) = {20, 22, 28, 26, 16, 14, 24, 18};
Volume(30) = {30};

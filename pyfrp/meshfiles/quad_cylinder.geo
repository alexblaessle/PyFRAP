volSize_px=22.0;
center_x=256;
center_y=256;
radius=294.103867936;
height=90.3332804037;


// Upper circle
Point(1) = {center_x, center_y, 0, volSize_px}; 
Point(2) = {center_x+radius, center_y, 0, volSize_px};
Point(3) = {center_x, center_y+radius, 0, volSize_px};

Point(4) = {center_x, center_y, -height, volSize_px}; 
Point(5) = {center_x+radius, center_y, -height, volSize_px};
Point(6) = {center_x, center_y+radius, -height, volSize_px};

//Circles describing upper circle
Circle(41) = {2, 1, 3}; 
Circle(42) = {5, 4, 6}; 

//Lines connecting circles 
Line(46) = {1,4};
Line(47) = {2,5};
Line(48) = {3,6};

Line(49) = {1,2};
Line(50) = {3,1};
Line(51) = {4,5};
Line(52) = {6,4};

Line Loop(54) = {49, 47, -51, -46};
Ruled Surface(54) = {54};
Line Loop(56) = {47, 42, -48, -41};
Ruled Surface(56) = {56};
Line Loop(58) = {46, -52, -48, 50};
Ruled Surface(58) = {58};
Line Loop(60) = {41, 50, 49};
Ruled Surface(60) = {60};
Line Loop(62) = {51, 42, 52};
Ruled Surface(62) = {62};
Surface Loop(64) = {56, 54, 60, 58, 62};
Volume(64) = {64};

volSize_px=23;
volSize_fine=11.5;
center_x=256.5;
center_y=256.5;
radius=110.863218314;
height=90.3332804037;
sl_refined=110.863218314;
sl_height=36.1333121615;
sl_width=3.9709623379;
volSize_fine=5.;

// Upper circle
Point(1) = {center_x, center_y, 0, volSize_px}; 
Point(2) = {center_x+radius, center_y, 0, volSize_px};
Point(3) = {center_x, center_y+radius, 0, volSize_px};

//Lower Circle
Point(4) = {center_x, center_y, -height, volSize_px}; 
Point(5) = {center_x+radius, center_y, -height, volSize_px};
Point(6) = {center_x, center_y+radius, -height, volSize_px};

//Slice Circle
Point(11) = {center_x, center_y, -sl_height, volSize_px}; 
Point(12) = {center_x+radius, center_y, -sl_height, volSize_px};
Point(13) = {center_x, radius+center_y, -sl_height, volSize_px};

//Circles describing upper circle
Circle(41) = {2, 1, 3}; 

//Circles describing lower circle
Circle(42) = {5, 4, 6}; 

//Circles describing slice circle
Circle(43) = {12, 11, 13}; 

//Lines connecting circles 
Line(46) = {1,4};
Line(47) = {2,5};
Line(48) = {3,6};

//Side lines of quadrants
Line(49) = {1,2};
Line(50) = {3,1};
Line(51) = {4,5};
Line(52) = {6,4};

//Side lines slice
Line(30) = {11,12};
Line(31) = {13,11};

//Slice Surface
Line Loop(35) = {43, 31, 30};
Plane Surface(36) = {35};

//Loops and Surfaces of quadrant
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

//Refined field on slice
Field[6] = Box;
Field[6].VIn = volSize_fine;
Field[6].VOut = volSize_px;
Field[6].XMin = center_x-sl_refined; 
Field[6].XMax = center_x+sl_refined;
Field[6].YMin = center_y-sl_refined;
Field[6].YMax = center_y+sl_refined;
Field[6].ZMin = -sl_height-sl_width;
Field[6].ZMax = -sl_height+sl_width;

Mesh.CharacteristicLengthExtendFromBoundary = 1;

//Set background field
Field[7] = Min;
Field[7].FieldsList = {6};
Background Field = 7;



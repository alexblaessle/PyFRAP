//Parameters
center_x=256.5;
center_y=256.5;
slice_height=-80.1333121615;
radius=250.563423869;
volSize_px=45;
volSize_fine=30;
sl_refined=110.863218314;
sl_height=36.1333121615;
sl_width=3.9709623379;


//formulas
outer_radius=(radius^2+slice_height^2)/(-2*slice_height);
inner_radius=outer_radius*1.1;
center_dist=Sqrt(inner_radius^2-outer_radius^2);
z_interc=-outer_radius;
x_interc=outer_radius;

//center inner ball
Point(1) = {center_x, center_y, -outer_radius-center_dist, volSize_px};

//top tip inner ball
Point(2) = {center_x, center_y, inner_radius-outer_radius-center_dist, volSize_px};

//center of intersection circle
Point(7) = {center_x, center_y, z_interc, volSize_px};

//4 Corner points of intersection circle
Point(8) = {center_x+x_interc, center_y, z_interc, volSize_px};
Point(9) = {center_x, center_y+x_interc, z_interc, volSize_px};
Point(10) = {center_x-x_interc, center_y, z_interc, volSize_px};
Point(11) = {center_x, center_y-x_interc, z_interc, volSize_px};

//Slice Circle
Point(3) = {center_x, center_y, slice_height, volSize_px}; 
Point(16) = {-radius+center_x, center_y, slice_height, volSize_px};
Point(17) = {center_x, center_y+radius, slice_height, volSize_px};
Point(18) = {center_x+radius, center_y, slice_height, volSize_px};
Point(19) = {center_x, -radius+center_y, slice_height, volSize_px};

//Circles describing slice circle
Circle(20) = {16, 3, 17}; 
Circle(21) = {17, 3, 18}; 
Circle(22) = {18, 3, 19}; 
Circle(23) = {19, 3, 16};

//Top tip outer circle
Point(29) = {center_x, center_y, 0, volSize_px};

//Circle of intersecting balls
Circle(12) = {8, 1, 2};
Circle(13) = {9, 1, 2};
Circle(14) = {10, 1, 2};
Circle(15) = {11, 1, 2};

//Arcs of balls
Circle(24) = {8, 7, 9};
Circle(25) = {9, 7, 10};
Circle(26) = {10, 7, 11};
Circle(27) = {11, 7, 8};
Circle(30) = {8, 7, 29};
Circle(31) = {9, 7, 29};
Circle(32) = {10, 7, 29};
Circle(33) = {11, 7, 29};

//Slice Surface
Line Loop(4) = {21, 22, 23, 20};
Plane Surface(5) = {4};

//Line Loops
Line Loop(35) = {33, -30, -27};
Ruled Surface(35) = {35};
Line Loop(37) = {32, -33, -26};
Ruled Surface(37) = {37};
Line Loop(39) = {31, -32, -25};
Ruled Surface(39) = {39};
Line Loop(41) = {24, 31, -30};
Ruled Surface(41) = {41};
Line Loop(43) = {12, -13, -24};
Ruled Surface(43) = {43};
Line Loop(45) = {13, -14, -25};
Ruled Surface(45) = {45};
Line Loop(47) = {14, -15, -26};
Ruled Surface(47) = {47};
Line Loop(49) = {15, -12, -27};
Ruled Surface(49) = {49};
Surface Loop(50) = {35, 37, 39, 41, 43, 49, 47, 45};
Volume(51) = {50};

//Box field
Field[6] = Box;
Field[6].VIn = volSize_fine;
Field[6].VOut = volSize_px;
Field[6].XMin = center_x-sl_refined; 
Field[6].XMax = center_x+sl_refined;
Field[6].YMin = center_y-sl_refined;
Field[6].YMax = center_y+sl_refined;
Field[6].ZMin = slice_height-sl_width;
Field[6].ZMax = slice_height+sl_width;

//Settings
Mesh.CharacteristicLengthExtendFromBoundary = 1;

//Background Field
Field[7] = Min;
Field[7].FieldsList = {6};
Background Field = 7;
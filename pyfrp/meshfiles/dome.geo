//Parameters
center_x=246.95;
center_y=234.28;
slice_height=-47.6515480548;
radius=243.61;
volSize_px=40.0;

//formulas
outer_radius=(radius^2+(slice_height)^2)/(-2*(slice_height));
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
//Top tip outer circle
Point(29) = {center_x, center_y, 0, volSize_px};
Circle(12) = {2, 1, 8};
Circle(13) = {2, 1, 9};
Circle(14) = {2, 1, 10};
Circle(15) = {2, 1, 11};
Circle(24) = {8, 7, 9};
Circle(25) = {9, 7, 10};
Circle(26) = {10, 7, 11};
Circle(27) = {11, 7, 8};
Circle(30) = {29, 7, 8};
Circle(31) = {29, 7, 9};
Circle(32) = {29, 7, 10};
Circle(33) = {29, 7, 11};
Line Loop(35) = {33, -30, 27};
Ruled Surface(35) = {35};
Line Loop(37) = {32, -33, 26};
Ruled Surface(37) = {37};
Line Loop(39) = {31, -32, 25};
Ruled Surface(39) = {39};
Line Loop(41) = {24, -31, 30};
Ruled Surface(41) = {41};
Line Loop(43) = {12, -13, 24};
Ruled Surface(43) = {43};
Line Loop(45) = {13, -14, 25};
Ruled Surface(45) = {45};
Line Loop(47) = {14, -15, 26};
Ruled Surface(47) = {47};
Line Loop(49) = {15, -12, 27};
Ruled Surface(49) = {49};
Surface Loop(50) = {35, 37, 39, 41, 43, 49, 47, 45};
Volume(51) = {50};

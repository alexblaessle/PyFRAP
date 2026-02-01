//Parameters
center_x=256;
center_y=256;
slice_height=-127.;
radius=305.25;
volSize_px=25.;

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

//circle points
Point(9) = {center_x, center_y+x_interc, z_interc, volSize_px};
Point(10) = {center_x+x_interc, center_y, z_interc, volSize_px};

//Top tip outer circle
Point(29) = {center_x, center_y, 0, volSize_px};

//Arcs of circles
Circle(13) = {9, 1, 2};
Circle(14) = {10, 1, 2};
Circle(25) = {9, 7, 10};
Circle(31) = {9, 7, 29};
Circle(32) = {10, 7, 29};

//Line connecting top of inner and outer ball
Line(51) = {2, 29};

//Surfaces
Line Loop(39) = {31, -32, -25};
Ruled Surface(39) = {39};
Line Loop(45) = {13, -14, -25};
Ruled Surface(45) = {45};
Line Loop(53) = {13, 51, -31};
Ruled Surface(53) = {53};
Line Loop(55) = {51, -32, 14};
Ruled Surface(55) = {55};

//Volume
Surface Loop(57) = {39, 53, 45, 55};
Volume(57) = {57};
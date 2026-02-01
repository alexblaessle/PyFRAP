//Parameters (for PyFRAP to work properly, the names of these three parameters must not be changed)
center_x=256;
center_y=256;
volSize_px=50;

//formulas (You can use all known mathematical expressions such as Sqrt, Sin,...)
some_parameter=Sqrt(center_x/center_y)+volSize_px

//Points
Point(id) = {x, y, z, volSize_px};

//Define Linear Edges here
Line(id) = {startPoint,endPoint};

//Define Circular Edges here
Circle(id) = {startPoint, centerPoint, endPoint};

//Define Line Loops here (note: endpoint and start points of successive edges must be matching)
Line Loop(id) = {ids of Edges};

//Define Surfaces here
Ruled Surface(id) = {id of Line Loop};

//Define Surface 
Surface Loop(id) = {ids, of, Surfaces};
Volume(id) = {id of Surface Loop};

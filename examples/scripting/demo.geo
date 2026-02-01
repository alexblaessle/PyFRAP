//vertices
Point(1)= {256.0,256.0,-15.0,10};
Point(2)= {556.0,256.0,-15.0,10};
Point(3)= {256.0,556.0,-15.0,10};
Point(4)= {-44.0,256.0,-15.0,10};
Point(5)= {256.0,-44.0,-15.0,10};

//lines

//arcs
Circle(1)= {2,1,3};
Circle(2)= {3,1,4};
Circle(3)= {4,1,5};
Circle(4)= {5,1,2};

//bSplines

//lineLoops
Line Loop(1)= {1,2,3,4};

//ruledSurfaces
Ruled Surface(1)= {1};

//surfaceLoops

//volumes

//fields


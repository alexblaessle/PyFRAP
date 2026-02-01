//vertices
Point(1)= {0,0,0,30.0};
Point(2)= {100,0,0,30.0};
Point(3)= {100,150,0,30.0};
Point(4)= {0,150,0,30.0};
Point(5)= {0,0,50,30.0};
Point(6)= {100,0,50,30.0};
Point(7)= {100,150,50,30.0};
Point(8)= {0,150,50,30.0};

//lines
Line(1)= {1,2};
Line(2)= {2,3};
Line(3)= {3,4};
Line(4)= {4,1};
Line(5)= {5,6};
Line(6)= {6,7};
Line(7)= {7,8};
Line(8)= {8,5};
Line(9)= {1,5};
Line(10)= {2,6};
Line(11)= {3,7};
Line(12)= {4,8};
Line(13)= {5,2};
Line(14)= {6,3};
Line(15)= {7,4};
Line(16)= {8,1};

//arcs

//bSplines

//lineLoops
Line Loop(1)= {5,-10,-13};
Line Loop(2)= {1,-13,-9};
Line Loop(3)= {6,-11,-14};
Line Loop(4)= {2,-14,-10};
Line Loop(5)= {7,-12,-15};
Line Loop(6)= {3,-15,-11};
Line Loop(7)= {8,-9,-16};
Line Loop(8)= {4,-16,-12};
Line Loop(9)= {8,5,6,7};
Line Loop(10)= {4,1,2,3};

//ruledSurfaces
Ruled Surface(1)= {1};
Ruled Surface(2)= {2};
Ruled Surface(3)= {3};
Ruled Surface(4)= {4};
Ruled Surface(5)= {5};
Ruled Surface(6)= {6};
Ruled Surface(7)= {7};
Ruled Surface(8)= {8};
Ruled Surface(9)= {9};
Ruled Surface(10)= {10};

//surfaceLoops

//volumes

//fields


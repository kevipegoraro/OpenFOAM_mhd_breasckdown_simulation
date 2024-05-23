// Gmsh project created on Thu Jun 18 00:29:38 2020
//Inputs
boxdim = 0.5;
refiner = 10;
gridsize = boxdim/refiner;

Point(1) = {0.425, 0, 0,gridsize};
Point(2) = {0.7, 0, 0,gridsize};
Point(3) = {0.825, 0.12, 0,gridsize};
Point(4) = {0.825, 0.36, 0,gridsize};
Point(5) = {0.7, 0.48, 0,gridsize};
Point(6) = {0.425, 0.48, 0,gridsize};  

Line(13) = {1,2};
Line(14) = {2,3};
Line(15) = {3,4};
Line(16) = {4,5};
Line(17) = {5,6};
Line(18) = {6,1};

Line Loop(25) = {13, 14, 15, 16, 17, 18};
Plane Surface(26) = 25;

Point(7) = {0.424882, 0, 0.009999077,gridsize};
Point(8) = {0.6998062, 0, 0.01647,gridsize};
Point(9) = {0.8247716, 0.12, 0.019412,gridsize};         
Point(10) = {0.8247716, 0.36, 0.019412,gridsize}; 
Point(11) = {0.6998062, 0.48, 0.01647,gridsize};
Point(12) = {0.424882, 0.48, 0.009999077,gridsize};

Line(19) = {8,7};
Line(20) = {9,8};
Line(21) = {10,9};
Line(22) = {11,10};
Line(23) = {12,11};
Line(24) = {7,12};

Line Loop(27) = {19, 20, 21, 22, 23, 24};
Plane Surface(28) = 27;

//made lines of edges
Line(29) = {1,7};  //fix
Line(30) = {8,2};  //fix
Line(31) = {3,9};  //fix
Line(32) = {10,4}; //fix
Line(33) = {5,11}; //fix
Line(34) = {12,6}; //fix

// made back plane
Line Loop(35) = {34, 18, 29, 24};
Plane Surface(36) = 35;
// made rigth_bottom plane
Line Loop(37) = {30, 14, 31, 20};
Plane Surface(38) = 37;
// made rigth_top plane
Line Loop(39) = {32, 16, 33, 22};
Plane Surface(40) = 39;
// made bottom plane
Curve Loop(40) = {29, -19, 30, -13};
Plane Surface(41) = {40};
// made top plane
Curve Loop(42) = {34, -17, 33, -23};
Plane Surface(43) = {42};
// made left plane
Curve Loop(44) = {31, -21, 32, -15};
Plane Surface(45) = {44};
// made the volume
Surface Loop(46) = {40, 45, 38, 41, 36, 43, 26, 28};
Volume(47) = {46};
// rename physical bouderies
Physical Surface("back") = {28};
Physical Surface("front") = {26};
Physical Surface("right") = {36};
Physical Surface("bottom") = {41};
Physical Surface("top") = {43};
Physical Surface("left") = {45};
Physical Surface("rigth_top") = {38};
Physical Surface("rigth_bottom") = {40};
Physical Volume("inter_mesh") = {47};

//Transfinite Line{5,6,7,8} = boxdim/gridsize+1;
//Transfinite Surface{10};
//Recombine Surface{14};



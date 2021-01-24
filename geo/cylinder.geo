r1 = 0.02;
r2 = 0.06;
lcr = 0.002;
nbElems_1 = r1/lcr+1;
nbElems_2 = r2/lcr+1;
nbElems_lateral = (r2-r1)/lcr+1;
// Inside
Point(1) = {0, 0, 0, lcr};
Point(2) = {0, r1, 0, lcr};
Point(3) = {r1, 0, 0, lcr};
Circle(1) = {2,1,3};
/* Transfinite Curve {1} = nbElems_1 Using Progression 1; */

// Outside
Point(4) = {0, r2, 0, lcr};
Point(5) = {r2, 0, 0, lcr};
Circle(2) = {5,1,4};

// Lateral
Line(3) = {3, 5};
Line(4) = {4,2};

/* Transfinite Curve {2} = nbElems_2 Using Progression 1; */
Transfinite Curve {3} = nbElems_lateral Using Progression 1;
Transfinite Curve {4} = nbElems_lateral Using Progression 1;
Line Loop(5) = {4, 1, 3, 2};

Plane Surface(1) = {5};

// Physical lines
Physical Line("bord_int") = {1};
Physical Line("bord_ext") = {2};
Physical Line("bord_inf") = {3};
Physical Line("bord_gauche") = {4};
Physical Surface("domaine") = {1};

// Mesh
Mesh 2;
//+

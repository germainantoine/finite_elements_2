Mesh.MshFileVersion = 2.2;
// definition du pas du maillage
h = 0.05;
// d finition des points (en 3D, raison pour laquelle il y a un 0 en z)
Point(1) = {1, 0, 0, h};
Point(2) = {6, 0, 0, h};
Point(3) = {6, 1, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {0, 0.5, 0, h};
Point(6) = {1, 0.5, 0, h};
// d finition des segments qui relient les points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
// d finition des contours ferm s
Line Loop(1) = {1,2,3,4,5,6};
// d finition des surfaces   partir contours ferm s
Plane Surface(1) = {1};
Physical Point(1) = {1,2,3,4,5,6};
Physical Line(1) = {4};
Physical Line(2) = {1,3,5,6};
Physical Line(3) = {2};
Physical Surface(1) = {1};

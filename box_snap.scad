// Box shaped snap connector

module box_snap(y, h, l, b, t=$fs, i_A=60, r_A=90){
// y = 3.28; // permissible deflection
// h = 6;    // thickness @ root
// l = 20;   // length of snap
// b = 10;   // width @ root

// t = 1;   // travel distance
// i_A = 45; // insert angle
// r_A = 60; // removal angle

i_t = y * tan(90 - i_A); // parallel travel on insertion 
r_t = y * tan(90 - r_A); // parallel travel on removal

rotate(a=90, 
v=[1,0,0])
linear_extrude(height=b, center=true)
polygon(points = [
    [0, 0],
    [r_t, y],
    [r_t + t, y],
    [r_t + t + i_t, 0],
    [r_t + t + i_t, -h],
    [-l, -h],
    [-l, 0]
    ]);
}

box_snap(y=3.28, h=6, l=50, b=10);
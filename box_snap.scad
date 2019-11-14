// Box shaped snap connector

// Tolerance
$fs = 0.2;

// Factor of safety
FOS = 10;

// Material properties
// ABS / PLA 
Sy = 35.0; // MPa
E = 2.3*1000; // MPa

function stress(
    F, // Force, N
    A  // Area, mm^2
    ) = F/A; // stress, MPa

function strain( 
    S, // Stress, MPa
    E, // Elastic modulus, MPa 
    ) = S / E ; // strain, %/100
 
// Permissible deflection 
function permissible_deflection(
    e_max, // maximum permissible strain, %/100
    l,     // length of arm, mm
    h,     // thickness @ root, mm
    K,     // geometric factor, 
    ) = (K*e_max)*pow(l,2)/h;   
function permissible_deflection_solve_for_h(e_max, l, y, K=0.67) =
    (K*e_max)*pow(l,2)/y; 
function permissible_deflection_solve_for_l(e_max, h, y, K=0.67) =
    pow(y*h/(K*e_max),0.5);

// Deflection force   
function deflection_force(
    e_max, // maximum permissible strain, %/100
    l, // length of arm, mm
    E, // elastic modulus, MPa
    Z // section modulus, mm^3
    ) = Z*E*e_max/l; // N
function deflection_force_solve_for_Z(e_max, l, E, P)
    = P*l/E/e_max;
function deflection_force_solve_for_l(e_max, E, Z, P)
    = P/E/e_max/Z;
    
// Mating force
function mating_force(
    P, // N, deflection force (perpendicular to insertion)
    A, // deg, mating angle
    mu, // Coefficent of friction
   ) =  ( friction_factor(A, mu) > 0 ? 
          friction_factor(A, mu) * P :
          1/0 ); // N
function mating_force_solve_for_A(P, mu, W)
    = friction_factor_solve_for_A(
        mu=mu, 
        K=W/P);
   
function friction_factor(A, mu)
    = (mu + tan(A))/(1 - mu*tan(A));
function friction_factor_solve_for_A(mu, K)
    = atan(K) - atan(mu);
 
// Moment of area
function moment_of_area_of_box(
    // Ix = integral(A=area){y^2} dA (about the x-axis)
    b, // mm, width @ root
    h, // mm, thickness @ root
    ) = b*pow(h,3)/12; // mm^4   
function moment_of_area_of_box_solve_for_b(h, I)
    = I/(pow(h,3)/12); // mm   
function moment_of_area_of_box_solve_for_h(b, I)
    = pow(12*(I/b),1/3); // mm

// Section modulus 
function section_modulus(
    c, // distance between neutral fibre and outer fibre
    I  // axial moment of inertia 
    // Z = I/c, 
    // where I = axial moment of inertia
    //       c = distance between neutral fibre and outer fibre
    ) = I/c; // mm^3
    
function section_modulus_of_box(
    b, // mm, width @ root 
    h, // mm, thickness @ root
    ) = section_modulus(
        c = h/2,
        I = moment_of_area_of_box(b=b, h=h) // mm^4
        ); 
function section_modulus_of_box_solve_for_b(h, Z)
    = moment_of_area_of_box_solve_for_b(
        h=h, 
        I=Z*(h/2)); // Z * c = I       
function section_modulus_of_box_solve_for_h(b, Z)
    // analytic solution to:
    // Z = bh^2/6
    = pow(Z/b*6, 0.5);

// Main
module snap_rectangle(
    // Main geometry
    y=false, // permissible deflection, mm
    h=false, // thickness @ root, mm
    l=false, // length of arm, mm
    b=false, // width @ root, mm
    // Head geometry
    t=$fs,  // travel distance, mm
    i_A=45, // insert angle, deg
    r_A=90, // removal angle, deg
    // Material properties
    Sy=Sy, // yield strength, MPa
    E=E,   // elastic modulus, MPa
    FOS=FOS, // Factor of safety, #
    // Deflection force
    P=false,  // Deflection force
    // Mating forces
    mu=0.5,   // coefficient of friction
    i_W=false, // insert mating force, N
    r_W=false, // removal mating force, N
    // Geometry 
    geometry=1, // 2d side geometry
    )
    {
      
    // K, geometric factor
    K = (geometry==1 ? 0.67 :
        (geometry==2 ? 1.09 :
        (geometry==3 ? 0.86 : 
          1/3)));  
    
    // strain (max permissible)
    e_max = strain(S=Sy/FOS, E=E);
    echo("Strain (max) (%/100) = ", e_max);
    
    // permissible deflection, mm   
    y = ( y ? y : permissible_deflection(
            e_max=e_max,
            l=l,
            h=h,
            K=K));
    echo("Permissible deflection (mm) = ", y);
     
    // thickness @ root, mm
    h = ( h ? h : permissible_deflection_solve_for_h(e_max, l, y, K));
    echo("Thickness @ root (mm) = ", h);
         
    // length of arm, mm
    l = ( l ? l : permissible_deflection_solve_for_l(e_max, h, y, K));
    echo("Length of arm (mm) = ", l);
    echo("Elongation (max) (mm) = ", l*e_max);
           
    // width @ root, mm
    b = ( b ? b : section_modulus_of_box_solve_for_b(
        h, 
        Z = P*l/E/e_max));
    echo("Width @ root (mm) = ", b);
    
    // deflection force
    Z = section_modulus_of_box(
        b, 
        h);
        
    P = ( P ? P : deflection_force(e_max, l, E, Z));
    echo("Deflection force (N) = ", P);
    
    // mating angles
    i_A = ( i_W ? mating_force_solve_for_A(P, mu, i_W) : i_A);
    r_A = ( r_W ? mating_force_solve_for_A(P, mu, r_W) : r_A);
    echo("Insertion angle (deg) = ", i_A);
    echo("Removal angle (deg) = ", r_A);
    
    // head insertion travel 
    i_t = y * tan(90 - i_A); // parallel travel on insertion 
    r_t = y * tan(90 - r_A); // parallel travel on removal
    
    // length of head, total
    echo("Length of head (mm) = ", t + i_t + r_t);
    echo("Total length (mm) = ", t + i_t + r_t + l);
    
    // mating force
    i_W = (i_W ? i_W : mating_force(P, i_A, mu));
    r_W = (r_W ? r_W : mating_force(P, r_A, mu));
    echo("Insertion force (N) = ", i_W);
    echo("Removal force (N) = ", r_W);
    
    // generate model
    if (geometry==1) {
        // Box snap   
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
    else if (geometry==2) {
        // Half snap h->h/2
        h2 = h/2;
        rotate(a=90, 
        v=[1,0,0])
        linear_extrude(height=b, center=true)
        polygon(points = [
            [0, 0],
            [r_t, y],
            [r_t + t, y],
            [r_t + t + i_t, 0],
            [r_t + t + i_t, -h2],
            [-l, -h],
            [-l, 0]
            ]);  
    }
    else if (geometry==3) {
        // Quarter snap b -> b/4
        b2 = b/4;
        points = [
            // Left
            [0, 0, b2/2], //0
            [r_t, y, b2/2],
            [r_t + t, y, b2/2],
            [r_t + t + i_t, 0, b2/2],
            [r_t + t + i_t, -h, b2/2],
            [0, -h, b2/2],
            [-l, -h, b/2],
            [-l, 0, b/2], //7
            // Right
            [0, 0, -b2/2], //8
            [r_t, y, -b2/2],
            [r_t + t, y, -b2/2],
            [r_t + t + i_t, 0, -b2/2],
            [r_t + t + i_t, -h, -b2/2],
            [0, -h, -b2/2],
            [-l, -h, -b/2],
            [-l, 0, -b/2] //15
            ];
        
        rotate(a=90, 
        v=[1,0,0])
        polyhedron(points=points, 
        faces=[[0,1,2,3,4,5,6,7], 
               [8,9,10,11,12,13,14,15],
               [0,1,9,8],
               [1,2,10,9],
               [2,3,11,10],
               [3,4,12,11],
               [4,5,13,12],
               [5,6,14,13],
               [6,7,15,14],
               [7,0,8,15]],
               convexity=10);
    }
    
    
    // Post-calcs   
    Area = (geometry==1 ? b*h :
            (geometry==2 ? b*h/2 :
            (geometry==3 ? b*h/4 : 
          b*h)));  
    msg1 = "Warning! Exceeds FOS";
    msg2 = "Warning!! Exceeds Yield Stress";

    echo("Insertion stress (MPa)", stress(i_W, Area));
    
    // assert( stress(i_W, Area) > Sy , msg2);
    // assert( stress(i_W, Area) > Sy/FOS, msg1);
    
    if ( stress(i_W, Area) > Sy ) {echo(msg2);}
    else if ( stress(i_W, Area) > Sy/FOS ) {echo(msg1);}
   
    
    echo("Removal stress (MPa)", stress(r_W, Area));
    if ( stress(r_W, Area) > Sy ) {echo(msg2);}
    else if ( stress(r_W, Area) > Sy/FOS ) {echo(msg1);}
    
        
}





// Demo
snap_rectangle(y=1, b=10, h=5, P=1, mu=0.5, geometry=1, t=1);

translate([0,-20,0])
snap_rectangle(y=1, b=10, h=5, P=1, mu=0.5, geometry=2, t=1);

translate([0,-40,0])
snap_rectangle(y=1, b=10, h=5, P=1, mu=0.5, geometry=3, t=1);
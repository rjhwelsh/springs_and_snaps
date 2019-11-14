// Box shaped snap connector

// Factor of safety
FOS = 10;

// material properties
// ABS / PLA 
Sy = 35.0; // MPa
E = 2.3*1000; // MPa

function strain( 
    S, // Stress, MPa
    E, // Elastic modulus, MPa 
    ) = S / E ; // strain, %/100
    
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
   
function section_modulus(
    c, // distance between neutral fibre and outer fibre
    I  // axial moment of inertia 
    // Z = I/c, 
    // where I = axial moment of inertia
    //       c = distance between neutral fibre and outer fibre
    ) = I/c; // mm^3
    
function moment_of_area_of_box(
    // Ix = integral(A=area){y^2} dA (about the x-axis)
    b, // mm, width @ root
    h, // mm, thickness @ root
    ) = b*pow(h,3)/12; // mm^4
    
function moment_of_area_of_box_solve_for_b(h, I)
    = I/(pow(h,3)/12); // mm
    
function moment_of_area_of_box_solve_for_h(b, I)
    = pow(12*(I/b),1/3); // mm

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
        
module box_snap(
    // Main geometry
    y=false, // permissible deflection, mm
    h=false, // thickness @ root, mm
    l=false, // length of arm, mm
    b=false, // width @ root, mm
    // Head geometry
    t=$fs,  // travel distance, mm
    i_A=60, // insert angle, deg
    r_A=90, // removal angle, deg
    // Material properties
    Sy=Sy, // yield strength, MPa
    E=E,   // elastic modulus, MPa
    FOS=FOS, // Factor of safety, #
    // Deflection force
    P=false,  // Deflection force
    // Geometry 
    K=0.67   // Geometric factor 
    )
    {
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
    
    // head insertion travel 
    i_t = y * tan(90 - i_A); // parallel travel on insertion 
    r_t = y * tan(90 - r_A); // parallel travel on removal
    
    // length of head, total
    echo("Length of head (mm) = ", t + i_t + r_t);
    echo("Total length (mm) = ", t + i_t + r_t + l);
    
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
    
    // generate model
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

box_snap(y=1, l=50, P=1);
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
    h     // thickness @ root, mm
    ) = (0.67*e_max)*pow(l,2)/h;
    
function permissible_deflection_solve_for_h(e_max, l, y) =
    (0.67*e_max)*pow(l,2)/y;
    
function permissible_deflection_solve_for_l(e_max, h, y) =
    pow(y*h/(0.67*e_max),0.5);

    
function deflection_force(
    e_max, // maximum permissible strain, %/100
    l, // length of arm, mm
    E, // elastic modulus, MPa
    Z // section modulus, mm^3
    ) = Z*E*e_max/l; // N

   
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

function section_modulus_of_box(
    b, // mm, width @ root 
    h, // mm, thickness @ root
    ) = section_modulus(
        c = h/2,
        I = moment_of_area_of_box(b=b, h=h) // mm^4
        ); 

        
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
    )
    {
    // strain (max permissible)
    e_max = strain(S=Sy/FOS, E=E);
    echo("Strain (max) (%/100) = ", e_max);
    
    // permissible deflection, mm   
    y = ( y ? y : permissible_deflection(
            e_max=e_max,
            l=l,
            h=h));
    echo("Permissible deflection (mm) = ", y);
     
    // width @ root, mm
    h = ( h ? h : permissible_deflection_solve_for_h(e_max, l, y));
    echo("Width @ root (mm) = ", h);
         
    // length of arm, mm
    l = ( l ? l : permissible_deflection_solve_for_l(e_max, h, y));
    echo("Length of arm (mm) = ", l);
    echo("Elongation (max) (mm) = ", l*e_max);  
    
    // head insertion travel 
    i_t = y * tan(90 - i_A); // parallel travel on insertion 
    r_t = y * tan(90 - r_A); // parallel travel on removal
    
    // length of head, total
    echo("Length of head (mm) = ", t + i_t + r_t);
    echo("Total length (mm) = ", t + i_t + r_t + l);
    
    // deflection force
    Z = section_modulus_of_box(
        b, 
        h);      
    P = ( P ? P : deflection_force(
                e_max,
                l,
                E, 
                Z));
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

box_snap(y=1, l=50, b=10);
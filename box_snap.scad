// Box shaped snap connector

// Factor of safety
FOS = 10;

// material properties
// ABS / PLA 
Sy = 35.0; // MPa
E = 2.3; // GPa

function strain( 
    S, // Stress, MPa
    E, // Elastic modulus, GPa 
    ) = S / ( E * 1000 ); // strain, %/100
    
function permissible_deflection(
    e_max, // maximum permissible strain, %/100
    l,     // length of arm, mm
    h     // thickness @ root, mm
    ) = (0.67*e_max)*pow(l,2)/h;
    
function permissible_deflection_solve_for_h(e_max, l, y) =
    (0.67*e_max)*pow(l,2)/y;
    
function permissible_deflection_solve_for_l(e_max, h, y) =
    pow(y*h/(0.67*e_max),0.5);

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
    E=E,   // elastic modulus, GPa
    FOS=FOS, // Factor of safety, #
    // Deflection force
    P=false,  // Deflection force
    Es=2000  // Secant modulus
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
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
function deflection_force_solve_for_E(e_max, l, Z, P)
    = P*l/e_max/Z;
    
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
function mating_force_solve_for_mu(P, A, W)
    = friction_factor_solve_for_mu(
        A=A,
        K=W/P);
   
function friction_factor(A, mu)
    = (mu + tan(A))/(1 - mu*tan(A));
function friction_factor_solve_for_A(mu, K)
    = atan(K) - atan(mu);
function friction_factor_solve_for_mu(A, K)
    = tan(max([atan(K) - A, 0]));
 
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
    e_max=false, // Maximum rated strain
    // Deflection force
    P=false,  // Deflection force
    // Mating forces
    mu=0.5,   // coefficient of friction
    i_W=false, // insert mating force, N
    r_W=false, // removal mating force, N
    // Geometry 
    geometry=1, // 2d side geometry
    // Metadata
    title="Untitled rectangular snap",
    border1="#############################",
    border2="_____________________________"
    )
    {
    echo(border1);
    echo(title);
    echo(border2);
        
    // K, geometric factor
    K = (geometry==1 ? 0.67 :
        (geometry==2 ? 1.09 :
        (geometry==3 ? 0.86 : 
          1/3)));  
    
    // strain (max permissible)
    e_max = (e_max ? e_max : strain(S=Sy/FOS, E=E));
    echo("Strain (max) (%/100) = ", e_max);
    FOS = ( e_max ? Sy/(E*e_max) : FOS );
    echo("Factor of Safety = ", FOS);
    echo(border2);
        
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
    echo(border2);
    
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
    echo(border2);


		module snap_head(h=h, a=0, b=b) {
				 // h = depth of snap head (from (0,0)) default:h
				 // a = angle of length i.e. atan(dh/L), 0=flat default:0
				 // b = width of snap head (across y-axis) default:b
				 let(dL = r_t + t + i_t,
						 points = [
									// Left
									[0, 0, b/2], //0
									[r_t, y, b/2],
									[r_t + t, y, b/2],
									[r_t + t + i_t, 0, b/2],
									[r_t + t + i_t, -h + dL*tan(a), b/2],
									[0, -h, b/2],
									// Right
									[0, 0, -b/2], //6
									[r_t, y, -b/2],
									[r_t + t, y, -b/2],
									[r_t + t + i_t, 0, -b/2],
									[r_t + t + i_t, -h + dL*tan(a), -b/2],
									[0, -h, -b/2]
									],
						 faces = [
									[0,1,7,6],
									[1,2,8,7],
									[2,3,9,8],
									[3,4,10,9],
									[4,5,11,10],

									[0,1,2,3,4,5],
									[6,7,8,9,10,11],

									[0,5,11,6]
									])
							rotate(a=90,
										 v=[1,0,0])
							polyhedron(points=points,
												 faces=faces,
												 convexity=10);
		}

		module snap_neck(h=h, l=l, b=b, h2=h, b2=b) {
				 // h = depth of snap neck (from (0,0)) default:h
				 // h2 = depth at snap head default:h
				 // l = length of snap neck default:l
				 // b = width of snap root (across y-axis) default:b
				 // b2 = width at snap head default:b
				 let(points = [// Left
									[0, 0, b2/2], //0
									[0, -h2, b2/2],
									[-l, -h, b/2],
									[-l, 0, b/2], //3
									// Right
									[0, 0, -b2/2], //4
									[0, -h2, -b2/2],
									[-l, -h, -b/2],
									[-l, 0, -b/2] //7
									],
						 faces=[
									[0,1,2,3],
									[4,5,6,7],
									[0,1,5,4],
									[1,2,6,5],
									[2,3,7,6],
									[3,0,4,7]]
							)
							rotate(a=90,
										 v=[1,0,0])
							polyhedron(points=points, faces=faces, convexity=10);
		}


		// generate model
		if (geometry==1) {
				 // Box snap
				 snap_neck();
				 snap_head();
		}
		else if (geometry==2) {
				 // Half snap h->h/2
				 h2=h/2;
				 snap_neck(h2=h2);
				 snap_head(h=h2, a=atan(h2/l));
		}
		else if (geometry==3) {
				 // Quarter snap b -> b/4
				 b2 = b/4;
				 snap_neck(b2=b2);
				 snap_head(h=h, b=b2);
		}


    // Post-calcs
    echo(border2);
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
    
    echo("Axial yield force (N) = ", Sy*Area);
    echo(border2);
}


// Demo
snap_rectangle(y=1, b=10, h=5, P=1, mu=0.5, geometry=1, t=1, title="Geometry 1");

translate([0,-20,0])
snap_rectangle(y=1, b=10, h=5, P=1, mu=0.5, geometry=2, t=1, title="Geometry 2");

translate([0,-40,0])
snap_rectangle(y=1, b=10, h=5, P=1, mu=0.5, geometry=3, t=1, title="Geometry 3");
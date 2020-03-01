// Box shaped snap connector

use<../extended_functions/extended_functions.scad>
// Provides sum(vList) for general vector summation

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

// Convenient modules
module rotate_extrude_around_x_axis(angle=360, convexity=2, $fn=$fn) {
		 rotate(a=-90, v=[0,0,1])
					rotate_extrude(angle=angle, convexity=convexity, $fn=$fn)
					rotate(a=90, v=[0,0,1])
					children();
}


// Main
module snap_rectangle(
		 // Main geometry
		 y=false, // permissible deflection, mm
		 h=false, // thickness @ root, mm
		 l=false, // total length of arm, mm
		 b=false, // width @ root, mm
		 // Head geometry
		 t=$fs,  // travel distance, mm
		 i_A=45, // insert angle, deg
		 r_A=90, // removal angle, deg
		 // Bend geometry
		 bend_r=false, // bend radius (to exterior of snap), mm
		 bend_l=false, // bend length, mm (does not include radius)
		 bend_angle=180, // bend angle, deg
		 bend_internal=false, // bend towards inside of snap (instead of outside)
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
		 echo("Total length (including snap head) (mm) = ", t + i_t + r_t + l);

		 // mating force
		 i_W = (i_W ? i_W : mating_force(P, i_A, mu));
		 r_W = (r_W ? r_W : mating_force(P, r_A, mu));
		 echo("Insertion force (N) = ", i_W);
		 echo("Removal force (N) = ", r_W);
		 echo(border2);

		 // Add 0 to front of bend_l, array only
		 bend_l = is_num(bend_l) ? bend_l : [ for (i=[-1:len(bend_l)-1]) i < 0 ? 0 : bend_l[i]];

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
									 [1,0,6,7],
									 [2,1,7,8],
									 [3,2,8,9],
									 [4,3,9,10],
									 [5,4,10,11],
									 [0,5,11,6],

									 [0,1,2,3,4,5],
									 [11,10,9,8,7,6]

									 ])
							 polyhedron(points=points,
													faces=faces,
													convexity=10);
		 }

		 module snap_neck(h=h, l=l, b=b, h2=h, b2=b, l_start=0, l_end=l) {
					// h = depth of snap neck (from (0,0)) default:h
					// l = length of snap neck default:l
					// b = width of snap root (across y-axis) default:b
					// Specify varying dimensions; h, b
					// h2 = depth at snap head default:h
					// b2 = width at snap head default:b
					// Specify a section of snap_neck only:
					// l_start = length at which to start default: 0
					// l_end   = length at which to end   default: l
					let(dh = h - h2,
							db = (b - b2)/2,

							ang_h = atan(dh/l),
							ang_b = atan(db/l),

							h_start = h2 + l_start*tan(ang_h),
							h_end = h2 + l_end*tan(ang_h),

							b_start = b2 + l_start*tan(ang_b),
							b_end = b2 + l_end*tan(ang_b),

							points = [// Left
									 [-l_start, 0, b_start/2], //0
									 [-l_start, -h_start, b_start/2],
									 [-l_end, -h_end, b_end/2],
									 [-l_end, 0, b_end/2], //3
									 // Right
									 [-l_start, 0, -b_start/2], //4
									 [-l_start, -h_start, -b_start/2],
									 [-l_end, -h_end, -b_end/2],
									 [-l_end, 0, -b_end/2] //7
									 ],
							faces=[
									 [0,1,2,3],
									 [7,6,5,4],
									 [4,5,1,0],
									 [2,1,5,6],
									 [3,2,6,7],
									 [0,3,7,4]]
							 )
							 polyhedron(points=points, faces=faces, convexity=10);
		 }

		 module snap_bend(r=bend_r, l=bend_l, a=bend_angle) {
					// Return radius bend at length, l.
					// r = radius of bend
					// l = length (at which bend begins)
					// a = angle of bend (180deg MAX)
					// children() ->
					// apply to snap_neck();
					translate([-l,0,0])
							 translate([0,-r,0])
							 rotate_extrude_around_x_axis(angle=a, convexity=10, $fn=$fn) /* TODO: Rotate snap_bend according to 'n' */
							 projection(cut=true)
							 translate([0,r,0])                //apply radius
							 rotate(a=90, v=[0,1,0])      //rotate around y-axis (stand-up on end)
							 translate([l,0,0])          //lower to set length
							 children();                  //snap_neck()
		 }

		 function bend_length(n=0, bend_l=bend_l, relative=false) =
					// Length along snap for each bend (excludes radius perimeters)
					// bend_l = length or list of lengths at which to bend snap connector
					// n = the nth bend of the snap connector
					// relative = return relative bend
					let( l_is_num = is_num(bend_l),
							 l_is_arr = is_list(bend_l)
							 )
					(l_is_num ?
					 ( relative ?
						 bend_l :
						 bend_l*n ) :
					 (l_is_arr ?
						let( bend_l_max = len(bend_l) - 1,
								 bend_l_last = bend_l[bend_l_max]
								 )
						// First bend length (always 0)
						( n == 0 ? 0 :
							// After exceeding last bend length value
							( n > bend_l_max ?
								( relative ? bend_length(n=n-bend_l_max, bend_l=bend_l_last, relative=true) :
									( bend_length(n=bend_l_max, bend_l=bend_l, relative=false) +
										bend_length(n=n-bend_l_max, bend_l=bend_l_last, relative=false))) :
								// Normal point in bend array
								( relative ?
									bend_l[n] :
									sum(slice(bend_l,0,n))
										 ))) :
						undef))
					;

		 function count_bends(bend_l, l) =
					// Return the number of bends specified
					// bend_l = length or list of length at which to bend snap connector
					// l = total length of the snap connector
					let( l_is_num = is_num(bend_l),
							 l_is_arr = is_list(bend_l)
							 )
					(l_is_num ? ceil(l/bend_l) :
					 (l_is_arr ?
						let( remainder = l - sum(bend_l),
								 bend_l_last = bend_l[len(bend_l)-1]
								 )
						len(bend_l) + count_bends(bend_l=bend_l_last, l=remainder) :
						undef))
					;

		 module snap_neck_translate_segment(r=bend_r, l=bend_l, a=bend_angle, n=0, h=h, dh=0, dr=0) {
					// Return length segment at from l_start to l_end
					// r = radius of bend
					// l = length (relative length of between each bend)
					// a = angle of bend
					// n = bend #

					// h = root thickness at snap_head()
					// dh = root thickness increase (per section 'n')

					// children() ->
					// apply to snap_neck(l_start, l_end);

					/* r0 = r; */
					/* r  = r0 + n*dr; */

					// Condition for reversing bend
					function reverse_condition(a, n) =
							 let(b =  (a*n)%360) b >= 180 ? true : false;

					// Angle restriction
					function bend_angle_restrict(a, n) =
							 let(b =  (a*n)%360) b < 180 ? b : 180 - b; // Maximum angle is 360, 0->180, 0->-180

					// Recursive translation functions
					// Translates length (based on length/angle; decoupled from bend translation)
					function relative_translation_for_length(l, a, n=0) =
							 (n==0 ? [0, 0, 0] :
								let(l_rel=bend_length(n=n, bend_l=l, relative=true),
										l_prev=bend_length(n=n-1, bend_l=l, relative=false),
										a_n = bend_angle_restrict(a, n-1),
										reverse = reverse_condition(a, n-1))
								(n>=1 ?
								 (reverse ?
									sum([ [l_rel*cos(a_n),
												 l_rel*sin(a_n),
												 0] ,
												relative_translation_for_length(l=l, a=a, n=n-1)]) :
									sum([ [-l_rel*cos(a_n),
												 -l_rel*sin(a_n),
												 0] ,
												relative_translation_for_length(l=l, a=a, n=n-1)])
											)
								 : undef ));

					// Centers radius section of bend at position
					function relative_translation_for_radius_to_center(r, a, n=0) =
							 (n==0 ? [0, 0, 0] :
								let(a_n = bend_angle_restrict(a, n))
								(n>=1 ? [r*sin(a_n),
												 -r*cos(a_n),
												 0] : undef ));

					// Translates bend (based on bend radius; decoupled from length translation)
					function relative_translation_for_bend(r, a, n=0, h=h, dh=dh, dr=dr) =
							 (n==0 ? [0, 0, 0] :
								let(a_n = bend_angle_restrict(a, n-1),
										reverse = reverse_condition(a, n-1),
										reverse_flip = reverse_condition(a, n) != reverse,

										// Adjust h_n according to thickness at each step
										h_n  = h + n*dh,
										dx_h = reverse_flip ? -h_n*sin(a) : 0,
										dy_h = reverse_flip ? h_n*cos(a) : 0,

										prev_n=n-1,

										hyp = 2*(r-dr)*sin(a/2)  // Hypotenuse for chord in bending arc
										 )

								(n>=1 ?

								 (reverse ?
									// Reverse orientation
									let(a_nr = abs(a_n), // Reverse angle (start from 0)
											dx  = -hyp*cos(a/2) - dx_h,  // dx is mirrored (negative)
											dy  = hyp*sin(a/2) + dy_h,
											dx_rot = dx*cos(-a_nr) - dy*sin(-a_nr),  // Rotate opposite direction (clockwise)
											dy_rot = dx*sin(-a_nr) + dy*cos(-a_nr))
									sum([ [-dx_rot, -dy_rot, 0] , // dx in reverse direction
												relative_translation_for_bend(r=r-dr, a=a, n=prev_n, dr=dr)]) :

									// Regular orientation
									let(dx  = hyp*cos(a/2) + dx_h,  // Relative change in position (relative coords)
											dy  = hyp*sin(a/2) + dy_h,
											dx_rot = dx*cos(a_n) - dy*sin(a_n),  // Rotate change in position (relative, aligned rotation)
											dy_rot = dx*sin(a_n) + dy*cos(a_n))
									sum([ [-dx_rot, -dy_rot, 0] ,
												relative_translation_for_bend(r=r-dr, a=a, n=prev_n, dr=dr)])
											)
								 : undef ));


					if (n==0) children(); // No bend
					else {
							 let(
										a_n = bend_angle_restrict(a, n),
										reverse = reverse_condition(a, n),
										tl = relative_translation_for_length(l=l, a=a, n=n),
										trc = relative_translation_for_radius_to_center(r=r, a=a, n=n),
										trb = relative_translation_for_bend(r=r, a=a, n=n),
										l_abs = bend_length(n=n, bend_l=l, relative=false)
										) {
										translate(tl)
												 translate(trc)
												 translate(trb)
												 if (reverse) {
															mirror([1,0,0])                 // Mirror snap_neck() and snap_bend()
																	 rotate(a=-a_n, v=[0,0,1])  // Rotate in opposite direction
																	 translate([0,r,0])
																	 translate([l_abs,0,0])
																	 children();
												 }
												 else {
															rotate(a=a_n, v=[0,0,1])
																	 translate([0,r,0])
																	 translate([l_abs,0,0])
																	 children();
												 }
							 }
					}
		 }

		 module snap_main(h2=h, b2=b) {

					module reverse_bend(bend_internal=bend_internal){
							 if (bend_internal) {
										translate([0,-h2,0])
												 mirror([0,1,0])
												 children(); }
							 else { children(); }
					}

					if (is_bool(bend_l) == false) {
							 /* echo("bend_l is a number!"); */
							 let(
										segments=count_bends(bend_l, l),
										colors=["red", "green", "blue", "purple", "orange"],
										dh = h - h2,
										dh_over_n = l/segments*tan(atan(dh/l)),
										bend_r_array = [ for (n = [0:segments]) bend_r ? bend_r : y/segments + h2 + (n+1)*dh_over_n]
										) {
										reverse_bend(bend_internal=bend_internal)
												 for (n = [0:segments-1]){
															let(
																	 bend_r0 = bend_r ? bend_r : y/segments + h2 + dh_over_n,
																	 bend_dr = bend_r ? 0 : dh_over_n,
																	 bend_ra = bend_r_array[n],  // adjusted bend radius

																	 // Length segments
																	 l_start = bend_length(n=n, bend_l=bend_l),
																	 l_end =  min(bend_length(n=n+1, bend_l=bend_l), l)
																	 )
																	 color(colors[n%len(colors)])
																	 snap_neck_translate_segment(r=bend_ra, l=bend_l, a=bend_angle, n=n,
																															 h=h2, dh=dh_over_n, dr=bend_dr) {
																	 snap_neck(l_start=l_start, l_end=l_end, h2=h2, b2=b2);
																	 if (n < segments - 1)  // Exclude radius on last segment
																				snap_bend(r=bend_ra, l=l_end, a=bend_angle) snap_neck(h2=h2, b2=b2);

																	 echo("l_start=", l_start);
																	 echo("l_end=", l_end);
															}
												 }
										// Bend information
										echo(border2);
										echo("Bend radius sequence, mm =", bend_r_array);
										echo("Bend radius, TOTAL=", sum(bend_r_array));
										echo(border2);
							 }}
					else {
							 snap_neck(l_start=0, l_end=l, h2=h2, b2=b2);
					}
					snap_head(h=h2, a=atan((h-h2)/l), b=b2);
		 }


		 // generate model
		 if (geometry==1) {
					// Box snap
					snap_main();
		 }
		 else if (geometry==2) {
					// Half snap h->h/2
					h2=h/2;
					snap_main(h2=h2);
		 }
		 else if (geometry==3) {
					// Quarter snap b -> b/4
					b2 = b/4;
					snap_main(b2=b2);
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
snap_rectangle(y=2, b=10, h=5, P=1, mu=0.5, geometry=1, t=1, title="Geometry 1", bend_l=[40,20,10], bend_angle=180);

translate([0,-100,0])
snap_rectangle(y=1, b=10, h=5, P=1, mu=0.5, geometry=2, t=1, title="Geometry 2", bend_l=20, bend_internal=true);

translate([0,-180,0])
snap_rectangle(y=1, b=10, h=5, P=1, mu=0.5, geometry=3, t=1, title="Geometry 3", bend_l=20, bend_angle=90);

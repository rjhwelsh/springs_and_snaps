// Material Mechanics
// Provides functions for material mechanics
// Copyright (C) 2020  rjhwelsh

// Default Material Properties
include<material/pla.scad>

// Mechanics
function stress(
		 F, // Force, N
		 A  // Area, m^2
		 ) = F/A; // stress, Pa

function strain(
		 S, // Stress, Pa
		 E=E, // Elastic modulus, Pa
		 ) = S / E ; // strain, %/100

function axial_deformation(
        F, // Force, N
        L, // Length, m
        A, // Area, m^2
        E=E // Elastic modulus, Pa
        ) = (F*L)/(A*E); // elongation

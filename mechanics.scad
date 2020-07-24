// Material Mechanics
// Provides functions for material mechanics
// Copyright (C) 2020  rjhwelsh

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

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

// TSUNER Copyright (C) 2025 Ernesto Guerrero Fernandez
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.



#pragma once

#include <iostream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define EPSILON      DBL_EPSILON
#define SGN(x)       ((fabs(x) < EPSILON) ? 0 : ((x) > 0) ? 1 : -1)
#define PI  	     3.14159265358979323846


// Sim parameters
#define G                9.81
#define EPS               1e-6 // Used for truncated velocity (divide)
#define VMAX         10.0 //max velocity
#define Re            6371e3  // Volumetric mean earth radius.


#define Manning     0 // Use Manning friction. Only use if SEDIMENT=0. NOT FINISHED, DO NOT USE.

#define UNITS       0.01  // Change units of the Initial condition to put them in SI (it doesn't change the bathymetry!).


#define BEDLOAD     1   // If 0, the bedload doesn't move. If 1, it moves normally.
#define SUSPENSION  1   // If 0, suspension load isn't included in the simulation. If 1, it is. If both suspension and bedload are 0, standard SW equations are reconvered.

#define order       2
#define RK          1  // If RK==1, perform Runge-Kutta. If 0, perfoms MUSCL-Hancook (second order ADER aka Cauchy-Kovaleskaya). Note that MUSCL-Hancook is not fully tested

#define WBC         1  // Wall type boundary conditions. 2=periodic, 1=wall, 0=free flow.
#define SPL         0  // Sponge layer. It has to be >=0. 0=no Sponge layer. Do not use if Read_files==0. Currently unavailable, DO NOT use

#define Read_files  1 // If 1, it reads the initial condition, bathymetry and boundary conditions from the files in data.dat. If 0, it reads the initial condition from the functions defined in cini.cxx.

// Likely, you do not want to modify parameters after this.

using namespace std;

// Size of the computational block for edges (Riemann solver)
#define NUM_THREADS_HORIZONTAL_EDGES 8
#define NUM_THREADS_VERTICAL_EDGES  8
// Size of the computational block for volumens
#define NUM_THREADS_HORIZONTAL_VOL 8
#define NUM_THREADS_VERTICAL_VOL  8

// Round up a/b to the bigger closest integrer
#define iDivUp(a,b)  (((a)%(b) != 0) ? ((a)/(b) + 1) : ((a)/(b)))

// Round up a/b to the smaller closest integrer
#define iDivDown(a,b)  ((a)/(b))

// Aligns a to the biggest inegrer cloeset to b and that its multiple of b
#define iAlignUp(a,b)  (((a)%(b) != 0) ? ((a)-(a) % (b)+(b)) : (a))

// Aligns a to the smallest inegrer cloeset to b and that its multiple of b
#define iAlignDown(a,b)  ((a)-(a) % (b))



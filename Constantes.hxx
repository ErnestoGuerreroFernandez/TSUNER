#ifndef _CONSTANTES_H_
#define _CONSTANTES_H_

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


// Par�metros
#define G                9.81
#define EPS               1e-6 // Used for truncated velocity (divide)
#define VMAX         10.0 //max velocity
#define Re            6371e3  // Volumetric mean earth radius, for Haversine formula.


#define model       0
#define solver      0

#define UNITS       0.01  // Change units of the Initial condition to put it in SI (it doesn't change the bathymetry!!).


#define order       2
#define RK          1

#define WBC         0  // Wall type boundary conditions. 2=periodic, 1=wall, 0=free flow.
#define SPL         0  // Sponge layer. It has to be >=0. 0=no Sponge layer.

#define Read_bathy  1  // If 0, it doesn't read the bathymetry input file and take whatever it's in cini as initial condition for the bathymetry. If 1, it reads the bathymetry data from the input file and ignore any geometry input data (such as xmin or npy).
#define Read_eta    1  // If 0, it reads the initial free surface and velocities from the cini file. If 1, it reads them from a file especified in data.
#define interpolate 1   // If 1, it will interpolate the input file data to respect the number of points npx and npy from the data file. Otherwise it will use one element per data.

using namespace std;

// Tama�o de un bloque en el procesamiento de las aristas
#define NUM_HEBRAS_ANCHO_ARI 8
#define NUM_HEBRAS_ALTO_ARI  8
// Tama�o de un bloque en el c�lculo del nuevo estado y el delta_T
#define NUM_HEBRAS_ANCHO_EST 8
#define NUM_HEBRAS_ALTO_EST  8

// Redondea a / b al mayor entero m�s cercano
#define iDivUp(a,b)  (((a)%(b) != 0) ? ((a)/(b) + 1) : ((a)/(b)))

// Redondea a / b al menor entero m�s cercano
#define iDivDown(a,b)  ((a)/(b))

// Alinea a al mayor entero m�ltiplo de b m�s cercano
#define iAlignUp(a,b)  (((a)%(b) != 0) ? ((a)-(a) % (b)+(b)) : (a))

// Alinea a al menor entero m�ltiplo de b m�s cercano
#define iAlignDown(a,b)  ((a)-(a) % (b))

#ifndef GPU_VAR
#define GPU_VAR
__device__ double d_OC;
__device__ double d_Opar;
__device__ double d_Q;
__device__ double d_DS;
__device__ double d_HQB;
__device__ double d_FR;
// __device__ int d_SPL;
#endif

#endif

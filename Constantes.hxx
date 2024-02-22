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


// Parametros
#define G                9.81
#define EPS               1e-6 // Used for truncated velocity (divide)
#define VMAX         10.0 //max velocity
#define Re            6371e3  // Volumetric mean earth radius.


#define Manning     0 // Use Manning friction. Only use if SEDIMENT=0. NOT FINISHED, DO NOT USE.

#define UNITS       0.01  // Change units of the Initial condition to put them in SI (it doesn't change the bathymetry!).


#define SEDIMENT    1   //If 0, the bedload doesn't move. The model is then equivalent to standar SW. If 1, it moves normally.
// #define bedload_formula 1

#define order       2
#define RK          1  // If RK==1, perform Runge-Kutta. If 0, perfoms MUSCL-Hancook (second order ADER aka Cauchy-Kovaleskaya)

#define WBC         0  // Wall type boundary conditions. 2=periodic, 1=wall, 0=free flow.
#define SPL         0  // Sponge layer. It has to be >=0. 0=no Sponge layer. Do not use if Read_files==0

#define Read_files  1  // If 1, it reads the initial condition, bathymetry and boundary conditions from the files in data.dat. If 0, it reads the initial condition from the functions defined in cini.cxx.

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

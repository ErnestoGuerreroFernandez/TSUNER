// GPU/GPU_kernels.cuh
#pragma once
#include "../TSUNER.h"
#include <cuda_runtime.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <netcdf.h>
#include <time.h>
#include "helper_timer.h"
#include "sharedmem.cuh"




// --- NetCDF Helper Functions ---
void saveNetCDF(double *V, double *vec, int num, int npx, int npy, int nvol, double current_time);
void write_time(int npx, int npy, int num, double sim_time);
void writeEtaNC(int npx, int npy, int num, double sim_time, double *vec);
void writevelNC(int npx, int npy, int num, int coord, double *vec);
void writeHNC(int npx, int npy, int num, double *vec);
void write_hmax(double *vec);
void write_eta_max(double *vec);
void write_bmin(double *vec);
void closeNC();
int initNC(char *input_prefix, double *bati, SimParams params);


// Volume Kernels
__global__ void GetDeltaT(double* d_V, double* d_deltaTV, SimParams params);
__global__ void states_updater(double* d_Vk, double* d_V0, double* d_V, int npx, int npy, int k, double epsilon_h, double invx, double invy, double *FL, double *FR, double *FU, double *FD, double *d, double current_time, double *d_CC, double *d_tCC, long unsigned int timeCC, double *dyrup, float *time_dyrup);
__global__ void hmax_calculator(double* d_hmax, double * d_V, double epsilon_h, int npx, int npy);
__global__ void bmin_calculator(double* d_bmin, double * d_V, int npx, int npy);
__global__ void eta_max_calculator(double* d_eta_max, double *d_V, double epsilon_h, int npx, int npy);
__global__ void set_inicializacion(double* d_hmax, int npx, int npy);
__global__ void set_inicializacion_bmin(double *d_bmin, int npx, int npy);
__global__ void set_inicializacion_eta_max(double *d_eta_max, int npx, int npy);



// Reconstructor
__global__ void reconstructor(double* d_V, double* reconsL, double* reconsR, double* reconsD, double* reconsU, double delta_T, double inv_dx, double inv_dy, SimParams params, double* d_Vk, double *d_V0);

// Volume Solver
__global__ void VolumeSolver (int *fixH, double current_time, SimParams params, double* reconsL, double* reconsR, double* reconsD, double* reconsU, double *FL, double *FR, double *FU, double *FD, double *d_CC, double *d_tCC, long unsigned int timeCC, double *tidesBd, float *tidesT, int tidesNPT);

// Edge Kernels
// __device__ void RiemannSolver_Vols(double hhl, double zb0, double2 qql, double hcl, double hhr, double zb1, double2 qqr, double hcr, int pos0, int pos1, const SimParams &params, double *FL, double *FR, int normal_x, int normal_y, int fixed_sediment);

// SedFlux
// __device__ void u_calculator(double h, double2 q, double epsilon_h, double2& u);
// __device__ void rh_eval(double h,double dh,double val, double& hr);
// __device__ void rzb_eval(double zc,double dz,double val, double& z);
// __device__ void rq_eval(double2 q,double2 dq, double val, double2& qr);
// __device__ double bedload(double h, double2 un, double db, const SimParams &params);
// __device__ void Cflux(double h, double2 u, double hc, double *F, int fixed_sediment, double db, const SimParams &params);
// __device__ void roe_average(double h0, double2 u0, double z0, double h1, double2 u1, double z1, double &hm, double2 &um, double &zm);
// __device__ void calc_eigen_aprox0(double h, double2 u, double *autovalores);
// __device__ void u_calculator_double(double h, double qx, double qy, double epsilon_h, double2& u);
// __device__ void calculate_erosion(double2 u, double modu, double Rp, double vs, const SimParams &params, double &Fe);
// __device__ void calculate_deposition(double c, double vs, const SimParams &params, double &Fd);
// __device__ void state_filter(double& h, double2& q, double &hc, double epsilon_h);
// __device__ void state_filter_double(double& h, double& qx, double& qy, double &hc, double epsilon_h);


template <typename T>
T MinReduction(T *d_in, int n);

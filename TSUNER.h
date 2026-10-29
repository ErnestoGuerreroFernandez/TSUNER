#pragma once

#include "Sim_parameters.hxx"
#include <sys/stat.h>
#include <fstream>
#include <netcdf.h>
// #include <vector_types.h>
#include <string.h>

struct SimParams {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double total_time, save_time;
    int npx, npy, nvol;
    double CFL;
    double Q, R;
    double porosity;
    double DS;
    double FR;
    double Opar;
    double OC;
    double Dsg;
    double HQB;
    double vmax;
    double epsilon_h;
    double dx;
    double dy;
};




#ifdef __cplusplus
extern "C" {
#endif
// Provides the double2 type to standard g++ without requiring full CUDA headers
int TSUNER(double *dVh, double2 *dVhu, double *dVH, int *fixH, SimParams params, char *input_prefix, double* d, double *CC, double *tCC, long unsigned int timeCC, double *dyrup, float *time_dyrup, double &tiempo, int CudaDevice, double *tides_bd, float *tides_time, int tidesNPT);

int checkCUDAsupport(int cudaDevice);
void initializeGPU();

#ifdef __cplusplus
}
#endif


// Problem Loader
int LoadProblemData(string input_file, string &input_prefix, double **dVh, double2 **dVhu, double **dVH, double **CC, double **tCC, int **fixH, double **dyrup, float **time_dyrup, unsigned long int &time_grid, SimParams &params, double **tides_bd, float **tides_time, int &tidesNPT);
bool FileExist(string fichero);
void freeMemoryCPU(double *dVh, double2 *dVhu, double* dVH, int *fixH, double *CC, double *tCC, double *dyrup, float *time_dyrup);


// Interpolations
void printProgress(double percentage);
void check_errGRD_inte(int iret);
void interpolate_CC(string file, double *CC,  double *tCC, int npx, int npy, double xmin, double xmax, double ymin, double ymax, double xmin_grid, double ymin_grid, long unsigned int time, double dx_grid, double dy_grid, double *H, int id, double SLR_correction);
void interpolate_spherics(double xmin, double xmax, double ymin, double ymax, double xmin_fun, double ymin_fun, int npx_fun, int npy_fun, double dx_fun, double dy_fun, int npx, int npy, double *fun, double *sol);
void fix_eta(double *p, double *H, int npx, int npy, double SLR_correction);
void spheric2cartesian(double xmin, double xmax, double ymin, double ymax, double *dx, double dy, int npx, int npy, double *fun, double *sol);
void spheric2cartesian_CC(double xmin, double xmax, double ymin, double ymax, double *dxE, double dyE, int npx, int npy, int time_grid, double *CC);
void dy_transform(double &val);
void dx_transform(double ymin, double dy, int tam, double xmax, double xmin, double *dx, double dx_cartesianas, double &dmax);


// Tides
void tides_prep(double *eta_aux, double **tides_bd, float **tides_time, int npx, int npy, string file_tides, double xmin, double xmax, double ymin, double ymax, int &npt);


// cini
double bathy(double x, double y);
double initial_h (double x, double y);
void velocity(double x , double y, double2& u);
int fixed_bed_CI(double x, double y);




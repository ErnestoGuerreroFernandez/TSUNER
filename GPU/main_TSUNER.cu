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



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include "helper_timer.h"
#include "Sim_parameters.hxx"
#include "Reduction_kernel_min.cu"
#include "Volume_kernel.cu"
#include "Reconstructor.cu"
#include "netcdf.cu"
#include "VolumeSolver.cu"



void saveNetCDF(double *V, double *vec, int num, int npx, int npy, int nvol, double current_time)
{
    double sim_time = current_time;
    double aux;
    double hn;
    int i, j, pos;
    
    int posv=0;
    
    write_time(npx, npy, num, sim_time);
    

    for (j=SPL; j<npy-SPL; j++) {
        for (i=SPL; i<npx-SPL; i++) {
            pos = (j*npx+i)*5;
            hn = V[pos] - V[pos + 3];
            vec[posv] = hn; 
            posv++;
        }
    }

    writeEtaNC(npx-2*SPL, npy-2*SPL, num, sim_time,vec);
    

    // u1x
    posv=0;
    for (j=SPL; j<npy-SPL; j++) {
        for (i=SPL; i<npx-SPL; i++){
            pos = (j*npx+i)*5;
            aux=(V[pos + 1])/(V[pos]+EPSILON);
            vec[posv] = aux;
            posv++;
        }
    }
    writevelNC(npx-2*SPL, npy-2*SPL, num,0, vec);
    
    // u1y
    posv=0;
    for (j=SPL; j<npy-SPL; j++) {
        for (i=SPL; i<npx-SPL; i++) {
            pos = (j*npx+i)*5;
            aux=(V[pos + 2])/(V[pos]+EPSILON);
            vec[posv] = aux;
            posv++;
        }
    }
    writevelNC(npx-2*SPL,npy-2*SPL,num,1,vec);
    
    // H
    posv=0;
    for (j=SPL; j<npy-SPL; j++) {
        for (i=SPL; i<npx-SPL; i++) {
            pos = (j*npx+i)*5;
            aux=(V[pos + 3]);
            vec[posv] = -aux;
            posv++;
        }
    }
    writeHNC(npx-2*SPL,npy-2*SPL,num,vec);
    
}

void freeMemoryGPU(double *d_V, double *d_V0, double* d_Vk, double* reconsL, double* reconsR, double *reconsD, double *reconsU, double* d_hmax, double *d_deltaTV, double *d_CC, double *d_tCC, double *FL, double *FRR, double *FU, double *FD, double *d_bmin, double *d_eta_max)
{
    cudaFree(d_hmax);
    cudaFree(d_bmin);
    cudaFree(d_eta_max);
    cudaFree(d_V);
    cudaFree(d_V0);
    cudaFree(d_Vk);
    
    cudaFree(d_deltaTV);
    
    cudaFree(reconsL);
    cudaFree(reconsR);
    cudaFree(reconsU);
    cudaFree(reconsD);
    
    cudaFree(FL);
    cudaFree(FRR);
    cudaFree(FU);
    cudaFree(FD);
    
    cudaFree(d_CC);
    cudaFree(d_tCC);
    
}

// Returns 0 if all good, 1 if not enough GPU memory and 2 if not enough CPU memory
extern "C" int TSUNER(double *dVh, double2 *dVhu, double *dVH, int *fixH, double xmin, double xmax, double ymin, double ymax, int npx, int npy, int nvol, double dx, double dy, char *input_prefix,  double total_time, double save_time, double CFL, double Q, double DS, double FRc, double Opar, double OC, double Dsg, double HQB, double vmax, double epsilon_h, double* d, double *CC, double *tCC, long unsigned int timeCC, double *dyrup, float *time_dyrup, double &sim_time, int CudaDevice, double *tides_bd, float *tides_time, int tidesNPT)
{
    
    double *d_V, *d_V0,*d_Vk;
    int *d_fixH;

    double *d_d;

    
    double *d_hmax, *d_bmin, *d_eta_max;
    double *reconsL, *reconsR, *reconsU, *reconsD;
    double *d_CC=nullptr;
    double *d_tCC=nullptr;
    double *d_dyrup;
    float *d_time_dyrup;
    double *d_tidesBd;
    float *d_tidesT;
    double *V;


    double *d_deltaTV;
    double *vec;
    int err;

    int niters = 0;
    
    double inv_dx = 1.0/dx;
    double inv_dy = 1.0/dy;
    
    int VolDoubleData_size  = nvol*sizeof(double);
    int VolDouble2Data_size = nvol*sizeof(double2);
    int sizeofdouble2 = sizeof(double2);
    int sizeofdouble = sizeof(double);
    V = new double [nvol*5];

    double current_time, delta_T;

    StopWatchInterface *timer = NULL;
    int iter;
    
    
    
    // Allocate GPU memory
    cudaSetDevice(CudaDevice);
    
    
    cudaMemcpyToSymbol(d_OC,&OC, sizeof(double));
    cudaMemcpyToSymbol(d_FR,&FRc, sizeof(double));
    cudaMemcpyToSymbol(d_Opar,&Opar, sizeof(double));
    cudaMemcpyToSymbol(d_DS,&DS, sizeof(double));
    cudaMemcpyToSymbol(d_Q,&Q, sizeof(double));
    cudaMemcpyToSymbol(d_HQB,&HQB, sizeof(double));
    cudaMemcpyToSymbol(d_Dsg,&Dsg, sizeof(double));
    
    
    err=cudaMalloc((void **)&d_hmax,VolDoubleData_size);
    err+=cudaMalloc((void **)&d_bmin,VolDoubleData_size);
    err+=cudaMalloc((void **)&d_eta_max, VolDoubleData_size);
    err+=cudaMalloc((void **)&d_V, 5*VolDoubleData_size);
    err+=cudaMalloc((void **)&d_V0, 5*VolDoubleData_size);
    err+=cudaMalloc((void **)&d_Vk, 5*VolDoubleData_size);

    cudaMalloc((void **)&d_d, SPL*sizeof(double));

    
    err+=cudaMalloc((void **)&d_fixH, nvol*sizeof(int));
    
    err+=cudaMalloc((void **)&d_deltaTV, VolDoubleData_size);
    
    err+=cudaMalloc((void **)&reconsL, 5*VolDoubleData_size);
    err+=cudaMalloc((void **)&reconsR, 5*VolDoubleData_size);
    err+=cudaMalloc((void **)&reconsU, 5*VolDoubleData_size);
    err+=cudaMalloc((void **)&reconsD, 5*VolDoubleData_size);
    
    double *FL, *FR, *FU, *FD;
    
    err += cudaMalloc((void **)&FL, VolDoubleData_size*5);
    err += cudaMalloc((void **)&FR, VolDoubleData_size*5);
    err += cudaMalloc((void **)&FU, VolDoubleData_size*5);
    err += cudaMalloc((void **)&FD, VolDoubleData_size*5);
    
    #if (Read_files==1)
    err+=cudaMalloc((void **)&d_CC, (2*npx+2*npy)*4*timeCC*sizeof(double));
    err+=cudaMalloc((void **)&d_tCC, timeCC*sizeof(double));
    #endif

    if(tidesNPT>0){
        err+=cudaMalloc((void **)&d_tidesBd, (2*npx+2*npy)*tidesNPT*sizeof(double));
        err+=cudaMalloc((void **)&d_tidesT, tidesNPT*sizeof(float));
    }

    err+=cudaMalloc((void **)&d_dyrup, nvol*40*sizeof(double));
    err+=cudaMalloc((void **)&d_time_dyrup, 40*sizeof(float));


    if (err >0) {
        cudaFree(d_hmax);
        cudaFree(d_bmin);
        cudaFree(d_eta_max);
        cudaFree(d_V);
        cudaFree(d_V0);
        cudaFree(d_Vk);
        cudaFree(d_fixH);
        
        cudaFree(d_deltaTV);
        
        cudaFree(reconsL);
        cudaFree(reconsR);
        cudaFree(reconsU);
        cudaFree(reconsD);
        
        cudaFree(FL);
        cudaFree(FR);
        cudaFree(FU);
        cudaFree(FD);
        
        cudaFree(d_CC);
        cudaFree(d_tCC);

        if (tidesNPT>0){
            cudaFree(d_tidesBd);
            cudaFree(d_tidesT);
        }

        cudaFree(d_d);

        return 1;
    }

    // Grid and block size for CUDA
    dim3 blockGridEst(iDivUp(npx, NUM_THREADS_HORIZONTAL_VOL), iDivUp(npy, NUM_THREADS_VERTICAL_VOL));
    dim3 blockGridEst_RiemannSolver(iDivUp(npx+1, NUM_THREADS_HORIZONTAL_VOL), iDivUp(npy+1, NUM_THREADS_VERTICAL_VOL));
    dim3 threadBlockEst(NUM_THREADS_HORIZONTAL_VOL, NUM_THREADS_VERTICAL_VOL);




    for (int i=0;i<nvol;i++){
        int pos = i*5;
        V[pos] = dVh[i];
        V[pos+1] = dVhu[i].x;
        V[pos+2] = dVhu[i].y;
        V[pos+3] = dVH[i];
        V[pos+4] = 0.0;
    }
    cudaMemcpy(d_V, V, 5*VolDoubleData_size, cudaMemcpyHostToDevice);



    // We copy data from CPU to GPU
    // cudaMemset(d_V,  0, 5*VolDoubleData_size);
    // cudaMemcpy(d_V, dVh, VolDoubleData_size, cudaMemcpyHostToDevice);
    // cudaMemcpy2D(d_V + nvol, sizeof(double), &(dVhu[0].x), sizeofdouble2, sizeof(double), nvol, cudaMemcpyHostToDevice);
    // cudaMemcpy2D(d_V + 2 * nvol, sizeof(double), &(dVhu[0].y), sizeofdouble2, sizeof(double), nvol, cudaMemcpyHostToDevice);
    // cudaMemcpy(d_V + 3*nvol, dVH, VolDoubleData_size, cudaMemcpyHostToDevice);

    cudaMemcpy(d_Vk, d_V, 5*VolDoubleData_size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(d_V0, d_V, 5*VolDoubleData_size, cudaMemcpyDeviceToDevice);

    cudaMemcpy(d_fixH, fixH, nvol*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_CC, CC, (2*npx+2*npy)*4*timeCC*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tCC, tCC, timeCC*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tidesBd, tides_bd, (2*npx+2*npy)*tidesNPT*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tidesT, tides_time, tidesNPT*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_d, d, SPL*sizeof(double), cudaMemcpyHostToDevice);
    

    vec = new double [nvol]();
    if (vec == NULL) {
        delete [] vec;
        
        freeMemoryGPU(d_V, d_V0, d_Vk, reconsL, reconsR, reconsD, reconsU, d_hmax,d_deltaTV, d_CC, d_tCC, FL, FR, FU, FD,d_bmin,d_eta_max);
        if (tidesNPT>0){
            cudaFree(d_tidesBd);
            cudaFree(d_tidesT);
        }
        return 2;
    }

    for (int i=0; i<nvol; i++) {
        vec[i] = dVH[i];
    }

    // cudaMemcpyToSymbol(d_OC,  &OC,   sizeof(double));
    // cudaMemcpyToSymbol(d_FR,  &FR,   sizeof(double));
    // cudaMemcpyToSymbol(d_Opar,&Opar, sizeof(double));
    // cudaMemcpyToSymbol(d_DS,  &DS,   sizeof(double));
    // cudaMemcpyToSymbol(d_Q,   &Q,    sizeof(double));
    // cudaMemcpyToSymbol(d_HQB, &HQB,  sizeof(double));
    // cudaMemcpyToSymbol(d_Dsg, &Dsg,  sizeof(double));
    double R, porosity;
    R = 1 + FRc/(Opar*G*DS);
    porosity = 1 - DS*sqrt(FRc/Opar)/Q;
    iter = initNC(input_prefix, npx, npy, xmin, xmax, ymin, ymax, total_time,vec, R, DS, FRc, porosity, OC);

    if (iter == 1) {
        freeMemoryGPU(d_V, d_V0, d_Vk, reconsL, reconsR, reconsD, reconsU, d_hmax,d_deltaTV, d_CC, d_tCC, FL, FR, FU, FD,d_bmin,d_eta_max);
        if (tidesNPT>0){
            cudaFree(d_tidesBd);
            cudaFree(d_tidesT);
        }
        return 2;
    }


    set_inicializacion<<<blockGridEst, threadBlockEst>>>(d_hmax,npx,npy);  // Set h_max to 0
    set_inicializacion_bmin<<<blockGridEst, threadBlockEst>>>(d_bmin,npx,npy);
    set_inicializacion_eta_max<<<blockGridEst, threadBlockEst>>>(d_eta_max, npx, npy);

    // CALCULO DEL DELTA_T INICIAL
    GetDeltaT<<<blockGridEst, threadBlockEst>>>(d_V, d_deltaTV, npx, npy, dx, dy, CFL, epsilon_h);
    delta_T = MinReduction<double>(d_deltaTV, nvol);
    fprintf(stdout, "Initial deltaT = %e sec\n", delta_T);


    iter = 1;
    sdkCreateTimer(&timer);
    sdkStartTimer(&timer);
    current_time = 0.0;
    double next_save_time = current_time;

    
    while (current_time < total_time) {
        if ((save_time >= 0.0) && (current_time >= next_save_time)) {
            cudaMemcpy(V, d_V, 5*VolDoubleData_size, cudaMemcpyDeviceToHost);
            cudaMemcpy(vec, d_hmax, VolDoubleData_size, cudaMemcpyDeviceToHost);
            write_hmax(vec);
            cudaMemcpy(vec, d_eta_max, VolDoubleData_size, cudaMemcpyDeviceToHost);
            write_eta_max(vec);
            saveNetCDF(V, vec, niters, npx, npy, nvol, current_time);
            niters++;
            next_save_time += save_time;
        }



        #if RK == 1
            for (int k=0; k<order; k++) {
        #else
            int k=0;
        {
        #endif

            // Reconstruction and frictions
            reconstructor<<<blockGridEst, threadBlockEst>>>(d_V, reconsL, reconsR, reconsD, reconsU, npx, npy, delta_T, inv_dx, inv_dy, epsilon_h, d_Vk, d_V0);

            // Solver
            VolumeSolver<<<blockGridEst_RiemannSolver, threadBlockEst>>>(d_fixH, npx, npy, current_time, epsilon_h, reconsL, reconsR, reconsD, reconsU, FL, FR, FU, FD, d_CC, d_tCC, timeCC, d_tidesBd, d_tidesT, tidesNPT);

            // Updater
            states_updater<<<blockGridEst, threadBlockEst>>>(d_Vk, d_V0, d_V, npx, npy, k, epsilon_h, delta_T*inv_dx, delta_T*inv_dy, FL, FR, FU, FD, d_d, current_time, d_CC, d_tCC, timeCC,d_dyrup,d_time_dyrup);
        }

        current_time += delta_T;

        hmax_calculator<<<blockGridEst, threadBlockEst>>>(d_hmax,d_V,epsilon_h,npx,npy);
        bmin_calculator<<<blockGridEst, threadBlockEst>>>(d_bmin,d_V,npx,npy);
        eta_max_calculator<<<blockGridEst, threadBlockEst>>>(d_eta_max, d_V, epsilon_h, npx, npy);
        GetDeltaT<<<blockGridEst, threadBlockEst>>>(d_V, d_deltaTV, npx, npy, dx, dy, CFL, epsilon_h);

        delta_T = MinReduction<double>(d_deltaTV, nvol);

        fprintf(stdout, "Iteration %3d, deltaT = %e sec, ", iter, delta_T);
        fprintf(stdout, "Time = %g sec\n", current_time);
        iter++;


    }

    cudaMemcpy(V, d_V, 5*VolDoubleData_size, cudaMemcpyDeviceToHost);
    saveNetCDF(V, vec, niters, npx, npy, nvol, current_time);

    cudaDeviceSynchronize();
    cudaMemcpy(vec, d_hmax, VolDoubleData_size, cudaMemcpyDeviceToHost);
    write_hmax(vec);
    cudaMemcpy(vec, d_eta_max, VolDoubleData_size, cudaMemcpyDeviceToHost);
    write_eta_max(vec);
    cudaMemcpy(vec, d_bmin, VolDoubleData_size, cudaMemcpyDeviceToHost);
    for (int i=0; i<nvol;i++) vec[i] = -vec[i];
    write_bmin(vec);
    sdkStopTimer(&timer);
    sim_time = sdkGetTimerValue(&timer)*0.001;  // In seconds

/*
    double maux=-1e30;
    for (int i=0;i<nvol;i++){
        // aux[i] -= dVH[i];
        // maux = max(aux[i],maux);
        printf("%f\n",vec[i]);
    }

    // cout<<maux<<endl;*/

    sdkDeleteTimer(&timer);
    freeMemoryGPU(d_V, d_V0, d_Vk, reconsL, reconsR, reconsD, reconsU, d_hmax,d_deltaTV, d_CC, d_tCC, FL, FR, FU, FD,d_bmin,d_eta_max);
    if (tidesNPT>0){
        cudaFree(d_tidesBd);
        cudaFree(d_tidesT);
    }
    delete [] vec;
    delete [] V;

    cudaFree(d_d);

    closeNC();

    return 0;
}



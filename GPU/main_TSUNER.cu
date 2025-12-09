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



void saveNetCDF(double *dVh, double2 *dVhu, double* dVH, double *vec, int num, int npx, int npy, int nvol, double current_time)
{
    double sim_time = current_time;
    double aux;
    double hn;
    int i, j, pos;
    
    int posv=0;
    
    write_time(npx, npy, num, sim_time);
    

    for (j=SPL; j<npy-SPL; j++) {
        pos = j*npx;
        for (i=SPL; i<npx-SPL; i++) {
            hn = dVh[pos + i] - dVH[pos + i];
            vec[posv] = hn; 
            posv++;
        }
    }

    writeEtaNC(npx-2*SPL, npy-2*SPL, num, sim_time,vec);
    

    // u1x
    posv=0;
    for (j=SPL; j<npy-SPL; j++) {
        pos = j*npx;
        for (i=SPL; i<npx-SPL; i++){
            aux=(dVhu[pos + i].x)/(dVh[pos+i]+EPSILON);
            vec[posv] = aux;
            posv++;
        }
    }
    writevelNC(npx-2*SPL, npy-2*SPL, num,0, vec);
    
    // u1y
    posv=0;
    for (j=SPL; j<npy-SPL; j++) {
        pos = j*npx;
        for (i=SPL; i<npx-SPL; i++) {
            aux=(dVhu[pos + i].y)/(dVh[pos+i]+EPSILON);
            vec[posv] = aux;
            posv++;
        }
    }
    writevelNC(npx-2*SPL,npy-2*SPL,num,1,vec);
    
    // H
    posv=0;
    for (j=SPL; j<npy-SPL; j++) {
        pos = j*npx;
        for (i=SPL; i<npx-SPL; i++) {
            aux=(dVH[pos + i]);
            vec[posv] = -aux;
            posv++;
        }
    }
    writeHNC(npx-2*SPL,npy-2*SPL,num,vec);
    
}

void freeMemoryGPU(double *d_dVh, double2 *d_dVhu, double* d_dVH, double *d_Ah, double *d_Azb, double2 *d_Ahu, double *d0_dVh, double2 *d0_dVhu, double *d0_dVH, double* reconsLh, double* reconsRh, double *reconsDh, double *reconsUh, double *reconsDzb, double *reconsUzb, double *reconsLzb, double *reconsRzb, double2* reconsDhu, double2* reconsUhu, double2 *reconsLhu, double2 *reconsRhu, double* d_hmax, double *d_deltaTV, double *d_CC, double *d_tCC, double *FL, double *FRR, double *FU, double *FD, double *d_bmin)
{
    cudaFree(d_hmax);
    cudaFree(d_bmin);
    cudaFree(d_dVh);
    cudaFree(d_dVhu);
    cudaFree(d_dVH);
    
    cudaFree(d_deltaTV);
    
    cudaFree(d0_dVh);
    cudaFree(d0_dVH);
    cudaFree(d0_dVhu);
    
    cudaFree(d_Ah);
    cudaFree(d_Azb);
    cudaFree(d_Ahu);
    
    cudaFree(reconsLh);
    cudaFree(reconsRh);
    cudaFree(reconsUh);
    cudaFree(reconsDh);
    cudaFree(reconsLzb);
    cudaFree(reconsRzb);
    cudaFree(reconsUzb);
    cudaFree(reconsDzb);
    cudaFree(reconsLhu);
    cudaFree(reconsRhu);
    cudaFree(reconsUhu);
    cudaFree(reconsDhu);
    
    cudaFree(FL);
    cudaFree(FRR);
    cudaFree(FU);
    cudaFree(FD);
    
    cudaFree(d_CC);
    cudaFree(d_tCC);
    
}

// Returns 0 if all good, 1 if not enough GPU memory and 2 if not enough CPU memory
extern "C" int TSUNER(double *dVh, double2 *dVhu, double *dVH, int *fixH, double xmin, double xmax, double ymin, double ymax, int npx, int npy, int nvol, double dx, double dy, char *input_prefix,  double total_time, double save_time, double CFL, double Q, double DS, double FR, double Opar, double OC, double HQB, double vmax, double epsilon_h, double* d, double *CC, double *tCC, long unsigned int timeCC, double *dyrup, float *time_dyrup, double &sim_time, int CudaDevice, double *tides_bd, float *tides_time, int tidesNPT)
{
    
    double2 *d_dVhu;
    double2 *d0_dVhu;
    double2 *d_Ahu;
    double *d_dVh,*d_dVH, *d0_dVh,*d_Ah,*d_Azb, *d0_dVH;
    int *d_fixH;

    double *d_d;

    
    
    double *d_hmax, *d_bmin;
    double *reconsLh, *reconsRh, *reconsUh, *reconsDh;
    double *reconsLzb, *reconsRzb, *reconsUzb, *reconsDzb;
    double2 *reconsLhu, *reconsRhu, *reconsUhu, *reconsDhu;
    double *d_CC;
    double *d_tCC;
    double *d_dyrup;
    float *d_time_dyrup;
    double *d_tidesBd;
    float *d_tidesT;

    double *d_deltaTV;
    double *vec;
    int err;

    int niters = 0;
    
    double inv_dx = 1.0/dx;
    double inv_dy = 1.0/dy;
    
    int VolDoubleData_size  = nvol*sizeof(double);
    int VolDouble2Data_size = nvol*sizeof(double2);
    
    double current_time, delta_T;

    StopWatchInterface *timer = NULL;
    int iter;
    
    
    
    // Allocate GPU memory
    cudaSetDevice(CudaDevice);
    
    
    cudaMemcpyToSymbol(d_OC,&OC, sizeof(double));
    cudaMemcpyToSymbol(d_FR,&FR, sizeof(double));
    cudaMemcpyToSymbol(d_Opar,&Opar, sizeof(double));
    cudaMemcpyToSymbol(d_DS,&DS, sizeof(double));
    cudaMemcpyToSymbol(d_Q,&Q, sizeof(double));
    cudaMemcpyToSymbol(d_HQB,&HQB, sizeof(double));
    
    
    err=cudaMalloc((void **)&d_hmax,VolDoubleData_size);
    err+=cudaMalloc((void **)&d_bmin,VolDoubleData_size);
    err+=cudaMalloc((void **)&d_dVh, VolDoubleData_size);
    err+=cudaMalloc((void **)&d_dVhu, VolDouble2Data_size);

    cudaMalloc((void **)&d_d, SPL*sizeof(double));

    
    err+=cudaMalloc((void **)&d_dVH, VolDoubleData_size);
    err+=cudaMalloc((void **)&d_fixH, nvol*sizeof(int));
    
    err+=cudaMalloc((void **)&d0_dVh, VolDoubleData_size);
    err+=cudaMalloc((void **)&d0_dVH, VolDoubleData_size);
    err+=cudaMalloc((void **)&d0_dVhu, VolDouble2Data_size);
    
    err+=cudaMalloc((void **)&d_deltaTV, VolDoubleData_size);
    
    err+=cudaMalloc( (void **)&d_Ah, VolDoubleData_size);
    err+=cudaMalloc( (void **)&d_Azb, VolDoubleData_size);
    err+=cudaMalloc( (void **)&d_Ahu, VolDouble2Data_size);
    

    err+=cudaMalloc((void **)&reconsLh, VolDoubleData_size);
    err+=cudaMalloc((void **)&reconsRh, VolDoubleData_size);
    err+=cudaMalloc((void **)&reconsUh, VolDoubleData_size);
    err+=cudaMalloc((void **)&reconsDh, VolDoubleData_size);
    
    err+=cudaMalloc((void **)&reconsLzb, VolDoubleData_size);
    err+=cudaMalloc((void **)&reconsRzb, VolDoubleData_size);
    err+=cudaMalloc((void **)&reconsUzb, VolDoubleData_size);
    err+=cudaMalloc((void **)&reconsDzb, VolDoubleData_size);
    
    err+=cudaMalloc((void **)&reconsLhu, VolDouble2Data_size);
    err+=cudaMalloc((void **)&reconsRhu, VolDouble2Data_size);
    err+=cudaMalloc((void **)&reconsDhu, VolDouble2Data_size);
    err+=cudaMalloc((void **)&reconsUhu, VolDouble2Data_size);
    
    double *FL, *FRR, *FU, *FD;
    
    err += cudaMalloc((void **)&FL, VolDoubleData_size*4);
    err += cudaMalloc((void **)&FRR, VolDoubleData_size*4);
    err += cudaMalloc((void **)&FU, VolDoubleData_size*4);
    err += cudaMalloc((void **)&FD, VolDoubleData_size*4);
    
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
        cudaFree(d_dVh);
        cudaFree(d_dVhu);
        cudaFree(d_dVH);
        cudaFree(d_fixH);
        
        cudaFree(d0_dVh);
        cudaFree(d0_dVH);
        cudaFree(d0_dVhu);
        
        cudaFree(d_deltaTV);
        
        cudaFree(d_Ah);
        cudaFree(d_Azb);
        cudaFree(d_Ahu);
        
        cudaFree(reconsLh);
        cudaFree(reconsRh);
        cudaFree(reconsUh);
        cudaFree(reconsDh);
        cudaFree(reconsLzb);
        cudaFree(reconsRzb);
        cudaFree(reconsUzb);
        cudaFree(reconsDzb);
        cudaFree(reconsLhu);
        cudaFree(reconsRhu);
        cudaFree(reconsUhu);
        cudaFree(reconsDhu);
        
        cudaFree(FL);
        cudaFree(FRR);
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
    dim3 threadBlockEst(NUM_THREADS_HORIZONTAL_VOL, NUM_THREADS_VERTICAL_VOL);

    // We copy data from CPU to GPU
    cudaMemcpy(d_dVh, dVh, VolDoubleData_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dVhu, dVhu, VolDouble2Data_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dVH, dVH, VolDoubleData_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_fixH, fixH, nvol*sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_CC, CC, (2*npx+2*npy)*4*timeCC*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tCC, tCC, timeCC*sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(d_tidesBd, tides_bd, (2*npx+2*npy)*tidesNPT*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tidesT, tides_time, tidesNPT*sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(d_d, d, SPL*sizeof(double), cudaMemcpyHostToDevice);
    
    cudaMemcpy(d0_dVh, d_dVh, VolDoubleData_size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(d0_dVH, d_dVH, VolDoubleData_size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(d0_dVhu, d_dVhu, VolDouble2Data_size, cudaMemcpyDeviceToDevice);
    
    cudaMemcpy(d_Ah, d_dVh, VolDoubleData_size,cudaMemcpyDeviceToDevice);
    cudaMemcpy(d_Azb, d_dVH, VolDoubleData_size,cudaMemcpyDeviceToDevice);
    cudaMemcpy(d_Ahu, d_dVhu, VolDouble2Data_size,cudaMemcpyDeviceToDevice);
    

    vec = new double [nvol]();
    if (vec == NULL) {
        delete [] vec;
        
        freeMemoryGPU(d_dVh, d_dVhu, d_dVH, d_Ah, d_Azb, d_Ahu, d0_dVh, d0_dVhu, d0_dVH, reconsLh, reconsRh, reconsDh, reconsUh, reconsDzb, reconsUzb, reconsLzb, reconsRzb, reconsDhu, reconsUhu, reconsLhu, reconsRhu,d_hmax,d_deltaTV, d_CC, d_tCC, FL, FRR, FU, FD,d_bmin);
        if (tidesNPT>0){
            cudaFree(d_tidesBd);
            cudaFree(d_tidesT);
        }
        return 2;
    }

    for (int i=0; i<nvol; i++) {
        vec[i] = dVH[i];
    }

    cudaMemcpyToSymbol(d_OC,  &OC,   sizeof(double));
    cudaMemcpyToSymbol(d_FR,  &FR,   sizeof(double));
    cudaMemcpyToSymbol(d_Opar,&Opar, sizeof(double));
    cudaMemcpyToSymbol(d_DS,  &DS,   sizeof(double));
    cudaMemcpyToSymbol(d_Q,   &Q,    sizeof(double));
    cudaMemcpyToSymbol(d_HQB, &HQB,  sizeof(double));
    double R, porosity;
    R = 1 + FR/(Opar*G*DS);
    porosity = 1 - DS*sqrt(FR/Opar)/Q;
    iter = initNC(input_prefix, npx, npy, xmin, xmax, ymin, ymax, total_time,vec, R, DS, FR, porosity, OC);

    if (iter == 1) {
        freeMemoryGPU(d_dVh, d_dVhu, d_dVH, d_Ah, d_Azb, d_Ahu, d0_dVh, d0_dVhu, d0_dVH, reconsLh, reconsRh, reconsDh, reconsUh, reconsDzb, reconsUzb, reconsLzb, reconsRzb, reconsDhu, reconsUhu, reconsLhu, reconsRhu,d_hmax,d_deltaTV, d_CC, d_tCC, FL, FRR, FU, FD,d_bmin);
        if (tidesNPT>0){
            cudaFree(d_tidesBd);
            cudaFree(d_tidesT);
        }
        return 2;
    }


    set_inicializacion<<<blockGridEst, threadBlockEst>>>(d_hmax,npx,npy);  // Set h_max to 0
    set_inicializacion_bmin<<<blockGridEst, threadBlockEst>>>(d_bmin,npx,npy);
    // CALCULO DEL DELTA_T INICIAL
    GetDeltaT<<<blockGridEst, threadBlockEst>>>(d_dVh, d_dVhu, d_dVH, d_deltaTV, npx, npy,dx, dy,CFL, epsilon_h);
    delta_T = MinReduction<double>(d_deltaTV, nvol);
    fprintf(stdout, "Initial deltaT = %e sec\n", delta_T);
    

    iter = 1;
    sdkCreateTimer(&timer);
    sdkStartTimer(&timer);
    current_time = 0.0;
    double next_save_time = current_time;

        double aux[npx*npy];
    for (int i=0;i<nvol;i++){
        aux[i] = dVH[i];
    }

    
    while (current_time < total_time) {
        if ((save_time >= 0.0) && (current_time >= next_save_time)) {
            cudaMemcpy(dVh, d_dVh, VolDoubleData_size, cudaMemcpyDeviceToHost);
            cudaMemcpy(dVH, d_dVH, VolDoubleData_size, cudaMemcpyDeviceToHost);
            cudaMemcpy(dVhu, d_dVhu, VolDouble2Data_size, cudaMemcpyDeviceToHost);
            cudaMemcpy(vec, d_hmax, VolDoubleData_size, cudaMemcpyDeviceToHost);
            write_hmax(vec);
            saveNetCDF(dVh, dVhu, dVH, vec, niters, npx, npy, nvol, current_time);
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
            reconstructor<<<blockGridEst, threadBlockEst>>>(d_dVh, d_dVhu, d_dVH, reconsLh, reconsRh, reconsDh, reconsUh, reconsLzb, reconsRzb, reconsDzb, reconsUzb, reconsLhu, reconsRhu, reconsDhu, reconsUhu, npx, npy, delta_T, inv_dx, inv_dy, epsilon_h,d_Ahu,d0_dVh,d0_dVhu);

            // Solver
            VolumeSolver(d_fixH, d_dVh, d_dVhu, d_dVH, d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, npx, npy, delta_T, dx, dy, current_time, epsilon_h, vmax, reconsLh, reconsRh, reconsDh, reconsUh, reconsLzb, reconsRzb, reconsDzb, reconsUzb, reconsLhu, reconsRhu, reconsDhu, reconsUhu, FL, FRR, FU, FD, d_CC, d_tCC, timeCC, d_tidesBd, d_tidesT, tidesNPT);

            // Updater
            states_updater<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, d_dVh,d_dVhu,d_dVH,npx, npy, k, epsilon_h, delta_T*inv_dx, delta_T*inv_dy, FL, FRR, FU, FD, d_d, current_time, d_CC, d_tCC, timeCC,d_dyrup,d_time_dyrup);
        }

        current_time += delta_T;

        hmax_calculator<<<blockGridEst, threadBlockEst>>>(d_hmax,d_dVh,npx,npy);
        bmin_calculator<<<blockGridEst, threadBlockEst>>>(d_bmin,d_dVH,npx,npy);
        GetDeltaT<<<blockGridEst, threadBlockEst>>>(d_dVh, d_dVhu, d_dVH, d_deltaTV, npx, npy,dx, dy,CFL, epsilon_h);

        delta_T = MinReduction<double>(d_deltaTV, nvol);

        fprintf(stdout, "Iteration %3d, deltaT = %e sec, ", iter, delta_T);
        fprintf(stdout, "Time = %g sec\n", current_time);
        iter++;


    }

    cudaMemcpy(dVh, d_dVh, VolDoubleData_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(dVH, d_dVH, VolDoubleData_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(dVhu, d_dVhu, VolDouble2Data_size, cudaMemcpyDeviceToHost);

    saveNetCDF(dVh, dVhu, dVH, vec, niters, npx, npy, nvol, current_time);

    cudaDeviceSynchronize();
    cudaMemcpy(vec, d_hmax, VolDoubleData_size, cudaMemcpyDeviceToHost);
    write_hmax(vec);
    cudaMemcpy(vec, d_bmin, VolDoubleData_size, cudaMemcpyDeviceToHost);
    for (int i=0; i<nvol;i++) vec[i] = -vec[i];
    write_bmin(vec);
    sdkStopTimer(&timer);
    sim_time = sdkGetTimerValue(&timer)*0.001;  // En segundos


    double maux=-1e30;
    for (int i=0;i<nvol;i++){
        aux[i] -= dVH[i];
        maux = max(aux[i],maux);
    }

    cout<<maux<<endl;

    sdkDeleteTimer(&timer);
    freeMemoryGPU(d_dVh, d_dVhu, d_dVH, d_Ah, d_Azb, d_Ahu, d0_dVh, d0_dVhu, d0_dVH, reconsLh, reconsRh, reconsDh, reconsUh, reconsDzb, reconsUzb, reconsLzb, reconsRzb, reconsDhu, reconsUhu, reconsLhu, reconsRhu,d_hmax,d_deltaTV, d_CC, d_tCC, FL, FRR, FU, FD,d_bmin);
    if (tidesNPT>0){
        cudaFree(d_tidesBd);
        cudaFree(d_tidesT);
    }
    delete [] vec;

    cudaFree(d_d);

    closeNC();

    return 0;
}



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



#include "GPU_kernels.cuh"


void freeMemoryGPU(double *d_V, double *d_V0, double* d_Vk, double* reconsL, double* reconsR, double *reconsD, double *reconsU, double* d_hmax, double *d_deltaTV, double *d_CC, double *d_tCC, double *FL, double *FR, double *FU, double *FD, double *d_bmin, double *d_eta_max, int *d_fixH, double* d_tidesBd, float* d_tidesT, int tidesNPT, double* d_d);

// Returns 0 if all good, 1 if not enough GPU memory and 2 if not enough CPU memory
extern "C" int TSUNER(double *dVh, double2 *dVhu, double *dVH, int *fixH, SimParams params, char *input_prefix, double* d, double *CC, double *tCC, long unsigned int timeCC, double *dyrup, float *time_dyrup, double &time_gpu, int CudaDevice, double *tides_bd, float *tides_time, int tidesNPT)
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

    double inv_dx = 1.0/params.dx;
    double inv_dy = 1.0/params.dy;

    int VolDoubleData_size  = params.nvol*sizeof(double);
    V = new double [params.nvol*5];

    double current_time, delta_T;

    StopWatchInterface *timer = NULL;
    int iter;
    
    
    // Allocate GPU memory
    cudaSetDevice(CudaDevice);
    
    err=cudaMalloc((void **)&d_hmax,VolDoubleData_size);
    err+=cudaMalloc((void **)&d_bmin,VolDoubleData_size);
    err+=cudaMalloc((void **)&d_eta_max, VolDoubleData_size);
    err+=cudaMalloc((void **)&d_V, 5*VolDoubleData_size);
    err+=cudaMalloc((void **)&d_V0, 5*VolDoubleData_size);
    err+=cudaMalloc((void **)&d_Vk, 5*VolDoubleData_size);

    err=+cudaMalloc((void **)&d_d, SPL*sizeof(double));

    
    err+=cudaMalloc((void **)&d_fixH, params.nvol*sizeof(int));
    
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
    err+=cudaMalloc((void **)&d_CC, (2*params.npx+2*params.npy)*4*timeCC*sizeof(double));
    err+=cudaMalloc((void **)&d_tCC, timeCC*sizeof(double));
    #endif

    if(tidesNPT>0){
        err+=cudaMalloc((void **)&d_tidesBd, (2*params.npx+2*params.npy)*tidesNPT*sizeof(double));
        err+=cudaMalloc((void **)&d_tidesT, tidesNPT*sizeof(float));
    }

    err+=cudaMalloc((void **)&d_dyrup, params.nvol*40*sizeof(double));
    err+=cudaMalloc((void **)&d_time_dyrup, 40*sizeof(float));


    if (err >0) {
        freeMemoryGPU(d_V, d_V0, d_Vk, reconsL, reconsR, reconsD, reconsU, d_hmax,d_deltaTV, d_CC, d_tCC, FL, FR, FU, FD,d_bmin,d_eta_max, d_fixH, d_tidesBd, d_tidesT, tidesNPT, d_d);
        return 1;
    }



    // Grid and block size for CUDA
    dim3 blockGridEst(iDivUp(params.npx, NUM_THREADS_HORIZONTAL_VOL), iDivUp(params.npy, NUM_THREADS_VERTICAL_VOL));
    dim3 blockGridEst_RiemannSolver(iDivUp(params.npx+1, NUM_THREADS_HORIZONTAL_VOL), iDivUp(params.npy+1, NUM_THREADS_VERTICAL_VOL));
    dim3 threadBlockEst(NUM_THREADS_HORIZONTAL_VOL, NUM_THREADS_VERTICAL_VOL);




    for (int i=0;i<params.nvol;i++){
        int pos = i*5;
        V[pos] = dVh[i];
        V[pos+1] = dVhu[i].x;
        V[pos+2] = dVhu[i].y;
        V[pos+3] = dVH[i];
        V[pos+4] = 0.0;
    }
    cudaMemcpy(d_V, V, 5*VolDoubleData_size, cudaMemcpyHostToDevice);



    cudaMemcpy(d_Vk, d_V, 5*VolDoubleData_size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(d_V0, d_V, 5*VolDoubleData_size, cudaMemcpyDeviceToDevice);

    cudaMemcpy(d_fixH, fixH, params.nvol*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_CC, CC, (2*params.npx+2*params.npy)*4*timeCC*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tCC, tCC, timeCC*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tidesBd, tides_bd, (2*params.npx+2*params.npy)*tidesNPT*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tidesT, tides_time, tidesNPT*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_d, d, SPL*sizeof(double), cudaMemcpyHostToDevice);
    

    vec = new double [params.nvol]();
    if (vec == NULL) {
        delete [] vec;
        freeMemoryGPU(d_V, d_V0, d_Vk, reconsL, reconsR, reconsD, reconsU, d_hmax,d_deltaTV, d_CC, d_tCC, FL, FR, FU, FD,d_bmin,d_eta_max, d_fixH, d_tidesBd, d_tidesT, tidesNPT, d_d);
        return 2;
    }

    for (int i=0; i<params.nvol; i++) {
        vec[i] = dVH[i];
    }


    iter = initNC(input_prefix, vec, params);

    if (iter == 1) {
        freeMemoryGPU(d_V, d_V0, d_Vk, reconsL, reconsR, reconsD, reconsU, d_hmax,d_deltaTV, d_CC, d_tCC, FL, FR, FU, FD,d_bmin,d_eta_max, d_fixH, d_tidesBd, d_tidesT, tidesNPT, d_d);
        return 2;
    }


    set_inicializacion<<<blockGridEst, threadBlockEst>>>(d_hmax,params.npx,params.npy);  // Set h_max to 0
    set_inicializacion_bmin<<<blockGridEst, threadBlockEst>>>(d_bmin,params.npx,params.npy);
    set_inicializacion_eta_max<<<blockGridEst, threadBlockEst>>>(d_eta_max, params.npx, params.npy);

    // Initial dt calculation
    GetDeltaT<<<blockGridEst, threadBlockEst>>>(d_V, d_deltaTV, params);
    delta_T = MinReduction<double>(d_deltaTV, params.nvol);
    fprintf(stdout, "Initial deltaT = %e sec\n", delta_T);


    iter = 1;
    sdkCreateTimer(&timer);
    sdkStartTimer(&timer);
    current_time = 0.0;
    double next_save_time = current_time;




    while (current_time < params.total_time) {
        if ((params.save_time >= 0.0) && (current_time >= next_save_time)) {
            cudaMemcpy(V, d_V, 5*VolDoubleData_size, cudaMemcpyDeviceToHost);
            cudaMemcpy(vec, d_hmax, VolDoubleData_size, cudaMemcpyDeviceToHost);
            write_hmax(vec);
            cudaMemcpy(vec, d_eta_max, VolDoubleData_size, cudaMemcpyDeviceToHost);
            write_eta_max(vec);
            saveNetCDF(V, vec, niters, params.npx, params.npy, params.nvol, current_time);
            niters++;
            next_save_time += params.save_time;
        }



        #if RK == 1
            for (int k=0; k<order; k++) {
        #else
            int k=0;
        {
        #endif

            // Reconstruction and frictions
            reconstructor<<<blockGridEst, threadBlockEst>>>(d_V, reconsL, reconsR, reconsD, reconsU, delta_T, inv_dx, inv_dy, params, d_Vk, d_V0);

            // Solver
            VolumeSolver<<<blockGridEst_RiemannSolver, threadBlockEst>>>(d_fixH, current_time, params, reconsL, reconsR, reconsD, reconsU, FL, FR, FU, FD, d_CC, d_tCC, timeCC, d_tidesBd, d_tidesT, tidesNPT);

            // Updater
            states_updater<<<blockGridEst, threadBlockEst>>>(d_Vk, d_V0, d_V, params.npx, params.npy, k, params.epsilon_h, delta_T*inv_dx, delta_T*inv_dy, FL, FR, FU, FD, d_d, current_time, d_CC, d_tCC, timeCC,d_dyrup,d_time_dyrup);
        }

        current_time += delta_T;

        hmax_calculator<<<blockGridEst, threadBlockEst>>>(d_hmax,d_V,params.epsilon_h,params.npx,params.npy);
        bmin_calculator<<<blockGridEst, threadBlockEst>>>(d_bmin,d_V,params.npx,params.npy);
        eta_max_calculator<<<blockGridEst, threadBlockEst>>>(d_eta_max, d_V, params.epsilon_h, params.npx, params.npy);
        GetDeltaT<<<blockGridEst, threadBlockEst>>>(d_V, d_deltaTV, params);

        delta_T = MinReduction<double>(d_deltaTV, params.nvol);

        fprintf(stdout, "Iteration %3d, deltaT = %e sec, ", iter, delta_T);
        fprintf(stdout, "Time = %g sec\n", current_time);
        iter++;


    }

    cudaMemcpy(V, d_V, 5*VolDoubleData_size, cudaMemcpyDeviceToHost);
    saveNetCDF(V, vec, niters, params.npx, params.npy, params.nvol, current_time);

    cudaDeviceSynchronize();
    cudaMemcpy(vec, d_hmax, VolDoubleData_size, cudaMemcpyDeviceToHost);
    write_hmax(vec);
    cudaMemcpy(vec, d_eta_max, VolDoubleData_size, cudaMemcpyDeviceToHost);
    write_eta_max(vec);
    cudaMemcpy(vec, d_bmin, VolDoubleData_size, cudaMemcpyDeviceToHost);
    for (int i=0; i<params.nvol;i++) vec[i] = -vec[i];
    write_bmin(vec);
    sdkStopTimer(&timer);
    time_gpu = sdkGetTimerValue(&timer)*0.001;  // In seconds



    sdkDeleteTimer(&timer);
    freeMemoryGPU(d_V, d_V0, d_Vk, reconsL, reconsR, reconsD, reconsU, d_hmax,d_deltaTV, d_CC, d_tCC, FL, FR, FU, FD,d_bmin,d_eta_max, d_fixH, d_tidesBd, d_tidesT, tidesNPT, d_d);
    delete [] vec;
    delete [] V;


    closeNC();


    return 0;
}


void freeMemoryGPU(double *d_V, double *d_V0, double* d_Vk, double* reconsL, double* reconsR, double *reconsD, double *reconsU, double* d_hmax, double *d_deltaTV, double *d_CC, double *d_tCC, double *FL, double *FR, double *FU, double *FD, double *d_bmin, double *d_eta_max, int *d_fixH, double* d_tidesBd, float* d_tidesT, int tidesNPT, double *d_d)
{
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
}




#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include "helper_timer.h"
#include "Constantes.hxx"
// #include "Arista_kernel.cu"
#include "Reduccion_kernel_min.cu"
#include "Volumen_kernel.cu"
#include "Reconstruye.cu"
#include "netcdf.cu"

#include "aristas_alternas.cu"
#include "aristas_volumenes.cu"

// __device__ double d_OC;
// __device__ double d_Opar;
// __device__ double d_Q;
// __device__ double d_DS;
// __device__ double d_HQB;
// __device__ double d_FR;
// __device__ int SPL;

void guardarNetCDF(double *dVh, double2 *dVhu, double* dVH, double *vec, int num, int npx, int npy, int nvol, double tiempo_act)
{
	double tiempo = tiempo_act;
	double aux;
        double hn;
	int i, j, pos;

        write_tiempo(npx, npy, num, tiempo);

	for (j=0; j<npy; j++) {
		pos = j*npx;
		for (i=0; i<npx; i++) {
			hn = dVh[pos + i] - dVH[pos + i];
			vec[j*npx + i] = (hn); 
                }
	}
    
	writeEtaNC(npx, npy, num, tiempo,vec);
    
    
    // u1x
    for (j=0; j<npy; j++) {
        pos = j*npx;
        for (i=0; i<npx; i++){
            aux=(dVhu[pos + i].x)/(dVh[pos+i]+EPSILON);
            vec[j*npx + i] = aux;
        }
    }
    writevelNC(npx, npy, num,0, vec);

    // u1y
    for (j=0; j<npy; j++) {
        pos = j*npx;
        for (i=0; i<npx; i++) {
            aux=(dVhu[pos + i].y)/(dVh[pos+i]+EPSILON);
            vec[j*npx + i] = aux;
        }
    }
    writevelNC(npx,npy,num,1,vec);
    
    // H
    for (j=0; j<npy; j++) {
        pos = j*npx;
        for (i=0; i<npx; i++) {
            aux=(dVH[pos + i]);
            vec[j*npx + i] = aux;
        }
    }
    writeHNC(npx,npy,num,vec);
	
}


void liberarMemoria(double *d_dVh, double2 *d_dVhu, double* d_dVH, double *d_Ah, double *d_Azb, double2 *d_Ahu, double *d0_dVh, double2 *d0_dVhu, double *d0_dVH, double* d_hx, double* d_hy, double2* d_hux, double2* d_huy, double* d_zbx, double* d_zby, double* d_hmax, double *d_deltaTV)
{
    cudaFree(d_hmax);
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
    
    cudaFree(d_hx);
    cudaFree(d_hy);
    cudaFree(d_hux);
    cudaFree(d_huy);
    cudaFree(d_zbx);
    cudaFree(d_zby);
    
}

// Devuelve 0 si todo ha ido bien, 1 si no hay memoria GPU suficiente, y 2 si no hay memoria CPU suficiente
extern "C" int Lava(double *dVh, double2 *dVhu, double *dVH, double xmin, double xmax, double ymin, double ymax, int npx, int npy, int nvol, double dx, double dy, char *prefijo,  double tiempo_tot, double tiempo_guardar, double CFL, double Q, double DS, double FR, double Opar, double OC, double HQB, double vmax, double epsilon_h, double* d, double hl, double Hl, double2 ql, double hr, double Hr, double2 qr, double ht, double Ht, double2 qt, double hb, double Hb, double2 qb, double &tiempo, int CudaDevice)
{
	int num_aristas_ver1, num_aristas_ver2;
	int num_aristas_hor1, num_aristas_hor2;
	
    double2 *d_dVhu;
    double2 *d0_dVhu;
	double2 *d_Ahu;
    double *d_dVh,*d_dVH, *d0_dVh,*d_Ah,*d_Azb, *d0_dVH; // Ah es una variable auxiliar para el RK
#if SPL>0
    double *d_d;
#endif

    
    double *d_hmax;
    double *d_hx, *d_hy, *d_zbx, *d_zby;
    double2 *d_hux, *d_huy;
    
    
    double *d_deltaTV;
	double *vec;
//  	cudaError_t err;
    int err;
	// Numero del estado que se va guardando
	int niters = 0;

	double inv_dx = 1.0/dx;
	double inv_dy = 1.0/dy;

	int tam_datosVolDouble  = nvol*sizeof(double);
	int tam_datosVolDouble2 = nvol*sizeof(double2);
    
	double tiempo_act, delta_T;
	double sig_tiempo_guardar = 0.0;
	StopWatchInterface *timer = NULL;
	int iter;
    

	// Numero de aristas verticales y horizontales
	num_aristas_ver1 = (int(npx/2) + 1)*npy;
	num_aristas_ver2 = (npx % 2 == 0) ? int(nvol/2) : num_aristas_ver1;
	num_aristas_hor1 = (int(npy/2) + 1)*npx;
	num_aristas_hor2 = (npy % 2 == 0) ? int(nvol/2) : num_aristas_hor1;
// 	num_aristas_ver1 = (int(npx/2) + 1)*npy;
// 	num_aristas_ver2 = (npx/2-1)*npy;
// 	num_aristas_hor1 = (int(npy/2) + 1)*npx;
// 	num_aristas_hor2 = (npy/2-1)*npx;
//     num_aristas_ver1 = npx/2*npy;
// 	num_aristas_ver2 = (npx/2)*npy;
// 	num_aristas_hor1 = npy/2*npx;
// 	num_aristas_hor2 = (npy/2)*npx;
	// Reservamos memoria en GPU
	cudaSetDevice(CudaDevice);

        
    cudaMemcpyToSymbol(d_OC,&OC, sizeof(double));
    cudaMemcpyToSymbol(d_FR,&FR, sizeof(double));
//     cudaMemcpyToSymbol(SPL,&SPL, sizeof(int));
    cudaMemcpyToSymbol(d_Opar,&Opar, sizeof(double));
    cudaMemcpyToSymbol(d_DS,&DS, sizeof(double));
    cudaMemcpyToSymbol(d_Q,&Q, sizeof(double));
    cudaMemcpyToSymbol(d_HQB,&HQB, sizeof(double));

    
    err=cudaMalloc((void **)&d_hmax,tam_datosVolDouble);
    err+=cudaMalloc((void **)&d_dVh, tam_datosVolDouble);
    err+=cudaMalloc((void **)&d_dVhu, tam_datosVolDouble2);
#if SPL>0
    cudaMalloc((void **)&d_d, SPL*sizeof(double));
#endif
    
    err+=cudaMalloc((void **)&d_dVH, tam_datosVolDouble);
    
    err+=cudaMalloc((void **)&d0_dVh, tam_datosVolDouble);
    err+=cudaMalloc((void **)&d0_dVH, tam_datosVolDouble);
    err+=cudaMalloc((void **)&d0_dVhu, tam_datosVolDouble2);
    
    err+=cudaMalloc((void **)&d_deltaTV, tam_datosVolDouble);
        
    err+=cudaMalloc( (void **)&d_Ah, tam_datosVolDouble);
    err+=cudaMalloc( (void **)&d_Azb, tam_datosVolDouble);
    err+=cudaMalloc( (void **)&d_Ahu, tam_datosVolDouble2);
    
    err+=cudaMalloc((void **)&d_hx, tam_datosVolDouble);
    err+=cudaMalloc((void **)&d_hy, tam_datosVolDouble);
    
	err+=cudaMalloc((void **)&d_hux, tam_datosVolDouble2);
    err+=cudaMalloc((void **)&d_huy, tam_datosVolDouble2);
	err+=cudaMalloc((void **)&d_zbx, tam_datosVolDouble);
    err+=cudaMalloc((void **)&d_zby, tam_datosVolDouble);
    
    double *FL, *FRR, *FU, *FD;
    
    err += cudaMalloc((void **)&FL, tam_datosVolDouble*4);
    err += cudaMalloc((void **)&FRR, tam_datosVolDouble*4);
    err += cudaMalloc((void **)&FU, tam_datosVolDouble*4);
    err += cudaMalloc((void **)&FD, tam_datosVolDouble*4);
    

    
	if (err >0) {
        cudaFree(d_hmax);
        cudaFree(d_dVh);
        cudaFree(d_dVhu);
        cudaFree(d_dVH);
        
		cudaFree(d0_dVh);
        cudaFree(d0_dVH);
        cudaFree(d0_dVhu);
        
        cudaFree(d_deltaTV);
        
		cudaFree(d_Ah);
        cudaFree(d_Azb);
        cudaFree(d_Ahu);

        cudaFree(d_hx);
		cudaFree(d_hy);
        cudaFree(d_hux);
        cudaFree(d_huy);
        cudaFree(d_zbx);
		cudaFree(d_zby);
#if SPL>0
        cudaFree(d_d);
#endif
		return 1;
	}

	// Tamano del grid y de bloque en el procesamiento de aristas
	dim3 blockGridVer1(iDivUp(num_aristas_ver1/npy, NUM_HEBRAS_ANCHO_ARI), iDivUp(npy, NUM_HEBRAS_ALTO_ARI));
	dim3 blockGridVer2(iDivUp(num_aristas_ver2/npy, NUM_HEBRAS_ANCHO_ARI), iDivUp(npy, NUM_HEBRAS_ALTO_ARI));
	dim3 blockGridHor1(iDivUp(npx, NUM_HEBRAS_ANCHO_ARI), iDivUp(num_aristas_hor1/npx, NUM_HEBRAS_ALTO_ARI));
	dim3 blockGridHor2(iDivUp(npx, NUM_HEBRAS_ANCHO_ARI), iDivUp(num_aristas_hor2/npx, NUM_HEBRAS_ALTO_ARI));
	dim3 threadBlock(NUM_HEBRAS_ANCHO_ARI, NUM_HEBRAS_ALTO_ARI);
	// Tamano del grid y bloques en el calculo del nuevo estado
	dim3 blockGridEst(iDivUp(npx, NUM_HEBRAS_ANCHO_EST), iDivUp(npy, NUM_HEBRAS_ALTO_EST));
	dim3 threadBlockEst(NUM_HEBRAS_ANCHO_EST, NUM_HEBRAS_ALTO_EST);
#if SPL>0
    // Tamano del grid y bloques en el calculo del sponge layer
	dim3 blockGridSPLBottomAndTop(iDivUp(npx, NUM_HEBRAS_ANCHO_EST), iDivUp(SPL, NUM_HEBRAS_ALTO_EST));
    dim3 blockGridSPLMid(iDivUp(npy-2*SPL, NUM_HEBRAS_ANCHO_EST), iDivUp(SPL, NUM_HEBRAS_ALTO_EST));
#endif
	// Copiamos los datos de CPU a GPU
	cudaMemcpy(d_dVh, dVh, tam_datosVolDouble, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dVhu, dVhu, tam_datosVolDouble2, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dVH, dVH, tam_datosVolDouble, cudaMemcpyHostToDevice);

#if SPL>0
	cudaMemcpy(d_d, d, SPL*sizeof(double), cudaMemcpyHostToDevice);
#endif
    
// 	cudaMemcpy(d0_dVh, d_dVh, tam_datosVolDouble, cudaMemcpyDeviceToDevice);
//     cudaMemcpy(d0_dVH, d_dVH, tam_datosVolDouble, cudaMemcpyDeviceToDevice);
//     cudaMemcpy(d0_dVhu, d_dVhu, tam_datosVolDouble2, cudaMemcpyDeviceToDevice);
	cudaMemcpy(d0_dVh, dVh, tam_datosVolDouble, cudaMemcpyHostToDevice);
    cudaMemcpy(d0_dVH, dVH, tam_datosVolDouble, cudaMemcpyHostToDevice);
    cudaMemcpy(d0_dVhu, dVhu, tam_datosVolDouble2, cudaMemcpyHostToDevice);
        
	
	// Inicializamos el acumulador
// 	cudaMemcpy(d_Ah, d_dVh, tam_datosVolDouble,cudaMemcpyDeviceToDevice);
//     cudaMemcpy(d_Azb, d_dVH, tam_datosVolDouble,cudaMemcpyDeviceToDevice);
// 	cudaMemcpy(d_Ahu, d_dVhu, tam_datosVolDouble2,cudaMemcpyDeviceToDevice);
    
    cudaMemcpy(d_Ah, dVh, tam_datosVolDouble, cudaMemcpyHostToDevice);
    cudaMemcpy(d_Azb, dVH, tam_datosVolDouble, cudaMemcpyHostToDevice);
	cudaMemcpy(d_Ahu, dVhu, tam_datosVolDouble2, cudaMemcpyHostToDevice);

    
    vec = new double [nvol]();
    if (vec == NULL) {
        delete [] vec;
        
        liberarMemoria(d_dVh, d_dVhu, d_dVH, d_Ah, d_Azb, d_Ahu, d0_dVh, d0_dVhu, d0_dVH, d_hx, d_hy, d_hux, d_huy, d_zbx, d_zby,d_hmax,d_deltaTV);
		return 2;
    }

    for (int i=0; i<nvol; i++) {
		vec[i] = dVH[i];
    }
    iter = initNC(prefijo, npx, npy, xmin, xmax, ymin, ymax, tiempo_tot,vec);
    if (iter == 1) {
        liberarMemoria(d_dVh, d_dVhu, d_dVH, d_Ah, d_Azb, d_Ahu, d0_dVh, d0_dVhu, d0_dVH, d_hx, d_hy, d_hux, d_huy, d_zbx, d_zby,d_hmax,d_deltaTV);
        return 2;
    }
    

//     set_inicializacion<<<blockGridEst, threadBlockEst>>>(flg,d_hmax,npx,npy);  // Pone a 0 flg y d_hmax
	// CALCULO DEL DELTA_T INICIAL
    obtenerDeltaTVolumenesGPU<<<blockGridEst, threadBlockEst>>>(d_dVh, d_dVhu, d_dVH, d_deltaTV, npx, npy,dx, dy,CFL, epsilon_h);
    delta_T = obtenerMinimoReduccion<double>(d_deltaTV, nvol);
	fprintf(stdout, "deltaT inicial = %e seg\n", delta_T);

	iter = 1;
	sdkCreateTimer(&timer);
	sdkStartTimer(&timer);
	tiempo_act = 0.0;
    
//     delta_T = 0.025463971798753 ;
    
    
//     for (int j=0; j<npy; j++) {
//         for (int i=0; i<npx; i++) {
//             int pos0 = (j-1)*npx + i;
//             int pos1 = j*npx + i;
//             if(j==0) pos0=(npy-1)*npx +i;
//             fprintf(stdout,"%.15f \n", (dVh)[pos0]-(dVh)[pos1]);
//             
//         }
//     }
   
    
	while (tiempo_act < tiempo_tot) {
		// Guardamos el estado actual, si procede
		// Inicio NetCDF
		if ((tiempo_guardar >= 0.0) && (tiempo_act >= sig_tiempo_guardar)) {
            

            cudaMemcpy(dVh, d_dVh, tam_datosVolDouble, cudaMemcpyDeviceToHost);
            cudaMemcpy(dVH, d_dVH, tam_datosVolDouble, cudaMemcpyDeviceToHost);
            cudaMemcpy(dVhu, d_dVhu, tam_datosVolDouble2, cudaMemcpyDeviceToHost);
            
            guardarNetCDF(dVh, dVhu, dVH, vec, niters, npx, npy, nvol, tiempo_act);
            niters++;

            sig_tiempo_guardar += tiempo_guardar;
        }

#if RK == 1
        for (int k=0; k<order; k++) {
#else
            int k=0;
            {
#endif

            reconstruye<<<blockGridEst, threadBlockEst>>>(d_dVh, d_dVhu, d_dVH, d_hx, d_hy, d_hux, d_huy, d_zbx,d_zby, npx, npy, delta_T, inv_dx, inv_dy, epsilon_h,d_Ahu);

//                 fprintf(stdout,"ANTES ---------------------------------------  \n");
//     borrar2<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, d_dVh,d_dVhu,d_dVH,npx, npy, k, epsilon_h, vmax,delta_T);   
//     cudaDeviceSynchronize();
//             aristas_alternas(d_dVh, d_dVhu, d_dVH, d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, npx, npy, delta_T, dx, dy, tiempo_act, epsilon_h, vmax, d_hx, d_hux, d_zbx, d_hy, d_huy, d_zby, num_aristas_hor1, num_aristas_hor2, num_aristas_ver1, num_aristas_ver2);
//     cudaDeviceSynchronize();
            
            aristas_volumenes(d_dVh, d_dVhu, d_dVH, d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, npx, npy, delta_T, dx, dy, tiempo_act, epsilon_h, vmax, d_hx, d_hux, d_zbx, d_hy, d_huy, d_zby, FL, FRR, FU, FD);
            
//     
//             fprintf(stdout,"DESPUES --------------------------------------- \n");
//     borrar2<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, d_dVh,d_dVhu,d_dVH,npx, npy, k, epsilon_h, vmax,delta_T);   
//     cudaDeviceSynchronize();
    

#if SPL>0            
            sponge_layerBottomAndTop<<<blockGridSPLBottomAndTop, threadBlock>>>(d_Ah, d_Ahu, d_Azb, npx, npy, d_d, hl, Hl, ql, hr, Hr, qr, ht, Ht, qt, hb, Hb, qb, tiempo_act);
            sponge_layerMid<<<blockGridSPLMid, threadBlock>>>(d_Ah, d_Ahu, d_Azb, npx, npy, d_d, hl, Hl, ql, hr, Hr, qr, tiempo_act);
#endif          
//             tiempo_act = 1000;
            actualiza_estados<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, d_dVh,d_dVhu,d_dVH,npx, npy, k, epsilon_h, vmax,delta_T);
//             cudaDeviceSynchronize();

        }

//         cudaMemcpy(d0_dVh,d_dVh, tam_datosVolDouble,cudaMemcpyDeviceToDevice);
//         cudaMemcpy(d0_dVH,d_dVH, tam_datosVolDouble,cudaMemcpyDeviceToDevice);
//         cudaMemcpy(d0_dVhu,d_dVhu,tam_datosVolDouble2,cudaMemcpyDeviceToDevice);
        
        // Actualizamos el tiempo actual
		tiempo_act += delta_T;

//         calculo_hmax<<<blockGridEst, threadBlockEst>>>(d_hmax,d_dVh,npx,npy);
        obtenerDeltaTVolumenesGPU<<<blockGridEst, threadBlockEst>>>(d_dVh, d_dVhu, d_dVH, d_deltaTV, npx, npy,dx, dy,CFL, epsilon_h);

        delta_T = obtenerMinimoReduccion<double>(d_deltaTV, nvol);

		fprintf(stdout, "Iteracion %3d, deltaT = %e seg, ", iter, delta_T);
		fprintf(stdout, "Tiempo = %g seg\n", tiempo_act);
		iter++;
	}
	
// 	fprintf(stdout,"ATENCION, SE CAMBIA EL DT INCIAL Y SE COMENTA EL CALCULO DEL DT INICIAL.");
	
    cudaMemcpy(dVh, d_dVh, tam_datosVolDouble, cudaMemcpyDeviceToHost);
    cudaMemcpy(dVH, d_dVH, tam_datosVolDouble, cudaMemcpyDeviceToHost);
    cudaMemcpy(dVhu, d_dVhu, tam_datosVolDouble2, cudaMemcpyDeviceToHost);
            
    guardarNetCDF(dVh, dVhu, dVH, vec, niters, npx, npy, nvol, tiempo_act);
        
	cudaDeviceSynchronize();
//         cudaMemcpy(vec, d_hmax, tam_datosVolDouble, cudaMemcpyDeviceToHost);
//         write_hmax(vec);
	sdkStopTimer(&timer);
	tiempo = sdkGetTimerValue(&timer)*0.001;  // En segundos

	sdkDeleteTimer(&timer);
    liberarMemoria(d_dVh, d_dVhu, d_dVH, d_Ah, d_Azb, d_Ahu, d0_dVh, d0_dVhu, d0_dVH, d_hx, d_hy, d_hux, d_huy, d_zbx, d_zby,d_hmax,d_deltaTV);
    
        delete [] vec;
#if SPL>0
        cudaFree(d_d);
#endif
        closeNC();

	return 0;
}


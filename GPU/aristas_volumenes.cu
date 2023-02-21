#ifndef _ARISTAS_VOLUMENES_H_
#define _ARISTAS_VOLUMENES_H_

#include <stdio.h>
#include "Arista_kernel.cu"
#include "Constantes.hxx"



__global__ void volumenes (double *d_dVh, double2 *d_dVhu, double *d_dVH, double *d_Ah, double2 *d_Ahu, double *d_Azb, double *d0_dVh, double2 *d0_dVhu, double *d0_dVH, int npx, int npy, double delta_T, double dx, double dy, double tiempo_act, double epsilon_h, double vmax, double *d_hx, double2 *d_hux, double *d_zbx, double *d_hy, double2 *d_huy, double *d_zby, double *FL, double *FR, double *FU, double *FD){

	int i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	int j = blockIdx.y*NUM_HEBRAS_ANCHO_EST  + threadIdx.y;
// 	int pi, pj, bp, bm;
    int pos0, pos1; 
    double h0, h1, zb0, zb1;
    double2 q0, q1;

    if ((i <= npx) && (j <= npy)) {
            if((i <= npx) && (j < npy)){
                // Riemann problems at vertical edges

                pos0 = i-1 + j*npx;
                pos1 = i + j*npx;
                
                
                if(i==0){
                    
#if WBC==2
                    
#else
                No periodicas no implementadas en volumenes
#endif
                }
                else{
                    if(i==npx){
                        pos1 = j*npx;
                        
                        h0 = d_dVh[pos0];
                        zb0 = d_dVH[pos0];
                        q0.x = d_dVhu[pos0].x;
                        q0.y = d_dVhu[pos0].y;
                        
                        h1 = d_dVh[pos1];
                        zb1 = d_dVH[pos1];
                        q1.x = d_dVhu[pos1].x;
                        q1.y = d_dVhu[pos1].y;
                        
                        pos0 *=4;
                        pos1*=4;
                        
                        procesarArista_Volumenes(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FL,FR,1,0);
                        
                    }
                    else{
                        h0 = d_dVh[pos0];
                        zb0 = d_dVH[pos0];
                        q0.x = d_dVhu[pos0].x;
                        q0.y = d_dVhu[pos0].y;
                        
                        h1 = d_dVh[pos1];
                        zb1 = d_dVH[pos1];
                        q1.x = d_dVhu[pos1].x;
                        q1.y = d_dVhu[pos1].y;
                        
                        pos0 *=4;
                        pos1*=4;
                        
                        procesarArista_Volumenes(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FL,FR,1,0);
                    }
                }
                
                
                
                
            }
            if((i < npx) && (j <= npy)){
                // Riemann problems at horizontal edges
                pos0 = i + (j-1)*npx;
                pos1 = i + j*npx;
                
                
                if(j==0){
                    
#if WBC==2
                    
#else
                No periodicas no implementadas en volumenes
#endif
                }
                else{
                    if(j==npy){
                        pos1 = i;
                        
                        h0 = d_dVh[pos0];
                        zb0 = d_dVH[pos0];
                        q0.x = d_dVhu[pos0].x;
                        q0.y = d_dVhu[pos0].y;
                        
                        h1 = d_dVh[pos1];
                        zb1 = d_dVH[pos1];
                        q1.x = d_dVhu[pos1].x;
                        q1.y = d_dVhu[pos1].y;
                        
                        pos0 *=4;
                        pos1*=4;
                        
                        procesarArista_Volumenes(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FD,FU,0,1);
                        
                    }
                    else{
                        h0 = d_dVh[pos0];
                        zb0 = d_dVH[pos0];
                        q0.x = d_dVhu[pos0].x;
                        q0.y = d_dVhu[pos0].y;
                        
                        h1 = d_dVh[pos1];
                        zb1 = d_dVH[pos1];
                        q1.x = d_dVhu[pos1].x;
                        q1.y = d_dVhu[pos1].y;
                        
                        pos0 *=4;
                        pos1*=4;
                        
                        procesarArista_Volumenes(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FD,FU,0,1);
                    }
                }
            }
        
    } 
    
    
    
    
}

__global__ void sum_contribuciones(double *d_Ah, double2 *d_Ahu, double *d_Azb, int npx, int npy, double delta_T, double dx, double dy, double epsilon_h, double vmax, double *FL, double *FR, double *FU, double *FD){
    
    int i,j;
    
    double invx, invy;
    
    invx = delta_T/dx;
    invy = delta_T/dy;

    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;
    
    
    if ((i<npx) && (j<npy)) {
        
        int pos = j*npx + i;
        int posF = pos*4;
        
        d_Ah[pos] -= invx*(FL[posF]+FR[posF])  + invy*(FD[posF]+FU[posF]);
        d_Ahu[pos].x -= invx*(FL[posF+1]+FR[posF+1])  + invy*(FD[posF+1]+FU[posF+1]);
        d_Ahu[pos].y -= invx*(FL[posF+2]+FR[posF+2]) + invy*(FD[posF+2]+FU[posF+2]);
        d_Azb[pos] -= invx*(FL[posF+3]+FR[posF+3])  + invy*(FD[posF+3]+FU[posF+3]);
        
        
        
    }
    
}

void aristas_volumenes(double *d_dVh, double2 *d_dVhu, double *d_dVH, double *d_Ah, double2 *d_Ahu, double *d_Azb, double *d0_dVh, double2 *d0_dVhu, double *d0_dVH, int npx, int npy, double delta_T, double dx, double dy, double tiempo_act, double epsilon_h, double vmax, double *d_hx, double2 *d_hux, double *d_zbx, double *d_hy, double2 *d_huy, double *d_zby, double *FL, double *FR, double *FU, double *FD)
{
    
    
    
    // Tamano del grid y de bloque en el procesamiento de aristas
//     dim3 threadBlock(NUM_HEBRAS_ANCHO_ARI, NUM_HEBRAS_ALTO_ARI);
    dim3 blockGridEst(iDivUp(npx, NUM_HEBRAS_ANCHO_EST), iDivUp(npy, NUM_HEBRAS_ALTO_EST));
    dim3 threadBlockEst(NUM_HEBRAS_ANCHO_EST, NUM_HEBRAS_ALTO_EST);
    
    
    volumenes<<<blockGridEst, threadBlockEst>>>(d_dVh, d_dVhu, d_dVH, d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, npx, npy, delta_T, dx, dy, tiempo_act, epsilon_h, vmax, d_hx, d_hux, d_zbx, d_hy, d_huy, d_zby, FL, FR, FU, FD);

    sum_contribuciones<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, npx, npy, delta_T, dx, dy, epsilon_h, vmax, FL, FR, FU, FD);
    
}
#endif 

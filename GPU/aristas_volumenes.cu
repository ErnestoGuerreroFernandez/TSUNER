#ifndef _ARISTAS_VOLUMENES_H_
#define _ARISTAS_VOLUMENES_H_

#include <stdio.h>
#include "Arista_kernel.cu"
#include "Constantes.hxx"



__device__ void interpolate_time(double *CC, double *tCC, long unsigned int timeCC, double tiempo_act, double &h1, double &zb1, double2 &q1, int npx, int npy, int pos){
    
    double dt;
    dt = tCC[1]-tCC[0];
    
    int t, tp;
    t = int((tiempo_act-tCC[0])/dt) ;
    t*=(t>0);
    tp = t+1;
    int posCC = t*(2*npx+2*npy)*4 + pos*4;
    int posCCp = tp*(2*npx+2*npy)*4 + pos*4;
    
    double m, b;
    m = (CC[posCCp]-CC[posCC])/dt;
    b = CC[posCC] - m*tCC[t];
    h1 = m*tiempo_act+b;
    h1*=UNITS;

    
    m = (CC[posCCp+3]-CC[posCC+3])/dt;
    b = CC[posCC+3] - m*tCC[t];
    zb1 = m*tiempo_act+b;
    
    h1 += zb1;    

    m = (CC[posCCp+1]-CC[posCC+1])/dt;
    b = CC[posCC+1] - m*tCC[t];
    q1.x = m*tiempo_act+b;    
    q1.x*= UNITS* h1   ;

    m = (CC[posCCp+2]-CC[posCC+2])/dt;
    b = CC[posCC+2] - m*tCC[t];
    q1.y = m*tiempo_act+b;   
    q1.y*= UNITS*h1    ;


    
}


__global__ void volumenes (double *d_dVh, double2 *d_dVhu, double *d_dVH, double *d_Ah, double2 *d_Ahu, double *d_Azb, double *d0_dVh, double2 *d0_dVhu, double *d0_dVH, int npx, int npy, double delta_T, double dx, double dy, double tiempo_act, double epsilon_h, double vmax, double* reconsLh, double* reconsRh, double* reconsDh, double* reconsUh, double* reconsLzb, double* reconsRzb, double* reconsDzb, double* reconsUzb, double2* reconsLhu, double2* reconsRhu, double2* reconsDhu, double2* reconsUhu, double *FL, double *FR, double *FU, double *FD, double *d_CC, double *d_tCC, long unsigned int timeCC){

	int i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	int j = blockIdx.y*NUM_HEBRAS_ANCHO_EST  + threadIdx.y;
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
                    pos0 = pos1;
                    pos1 = -1;
                        
                    h1 = h0 = reconsLh[pos0];
                    zb1 = zb0 = reconsLzb[pos0];
                    q1.x = q0.x = reconsLhu[pos0].x;
                    q1.y = q0.y = reconsLhu[pos0].y;
                        
#if WBC==1
                    q1.x = -q0.x;
#endif
#if (SPL<=0 && Read_files==1)  // Enforcing
                    if(tiempo_act<d_tCC[timeCC-2]){
                        interpolate_time(d_CC,d_tCC,timeCC,tiempo_act,h1,zb1,q1,npx,npy,2*npx+npy+j);
                    }
                    else{ // Salida libre
                        q1.x = q0.x;
                    }
#endif
                    pos0 *=4;

                    procesarArista_Volumenes(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FL,FR,-1,0);                  
#endif
                }
                else{
                    if(i==npx){
#if WBC==2
                        pos1 = j*npx;
                        h0 = reconsRh[pos0];
                        zb0 = reconsRzb[pos0];
                        q0.x = reconsRhu[pos0].x;
                        q0.y = reconsRhu[pos0].y;
                        
                        h1 = reconsLh[pos1];
                        zb1 = reconsLzb[pos1];
                        q1.x = reconsLhu[pos1].x;
                        q1.y = reconsLhu[pos1].y;
                        
                        pos1 *=4;
                        
#else
                        pos1 = -2;
                        
                        h1 = h0 = reconsRh[pos0];
                        zb1 = zb0 = reconsRzb[pos0];
                        q1.x = q0.x = reconsRhu[pos0].x;
                        q1.y = q0.y = reconsRhu[pos0].y;
                        
#if WBC==1
                        q1.x = -q0.x;
#endif
                        
#endif
#if (SPL<=0 && Read_files==1)  // Enforcing
                    if(tiempo_act<d_tCC[timeCC-2]){
                        interpolate_time(d_CC,d_tCC,timeCC,tiempo_act,h1,zb1,q1,npx,npy,npx+j);
                    }
                    else{ // Salida libre
                        q1.x = q0.x;
                    }
#endif

                        pos0 *=4;
                        procesarArista_Volumenes(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FL,FR,1,0);
                        
                    }
                    else{
                        h0 = reconsRh[pos0];
                        zb0 = reconsRzb[pos0];
                        q0.x = reconsRhu[pos0].x;
                        q0.y = reconsRhu[pos0].y;
                        
                        h1 = reconsLh[pos1];
                        zb1 = reconsLzb[pos1];
                        q1.x = reconsLhu[pos1].x;
                        q1.y = reconsLhu[pos1].y;
   
                        pos0 *=4;
                        pos1 *=4;
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
                    pos0 = pos1;
                    pos1 = -1;
                    h1 = h0 = reconsDh[pos0];
                    zb1 = zb0 = reconsDzb[pos0];
                    q1.x = q0.x = reconsDhu[pos0].x;
                    q1.y = q0.y = reconsDhu[pos0].y;
#if WBC==1
                    q1.y = -q0.y;
#endif
                    
#if (SPL<=0 && Read_files==1)  // Enforcing
                    if(tiempo_act<d_tCC[timeCC-2]){
                        interpolate_time(d_CC,d_tCC,timeCC,tiempo_act,h1,zb1,q1,npx,npy,i);
                    }
                    else{ // Salida libre
                        q1.y = q0.y;
                    }
#endif

                    
                    pos0 *=4;
                    procesarArista_Volumenes(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FD,FU,0,-1);
#endif
                }
                else{
                    if(j==npy){
#if WBC==2                        
                        pos1 = i;
                        
                        h0 = reconsUh[pos0];
                        zb0 = reconsUzb[pos0];
                        q0.x = reconsUhu[pos0].x;
                        q0.y = reconsUhu[pos0].y;
                        
                        h1 = reconsDh[pos1];
                        zb1 = reconsDzb[pos1];
                        q1.x = reconsDhu[pos1].x;
                        q1.y = reconsDhu[pos1].y;
                        
                        pos1*=4;
                        
#else
                        pos1 = -2;
                        h1 = h0 = reconsUh[pos0];
                        zb1 = zb0 = reconsUzb[pos0];
                        q1.x = q0.x = reconsUhu[pos0].x;
                        q1.y = q0.y = reconsUhu[pos0].y;
#if WBC==1
                        q1.y = -q0.y;
#endif
                        
#endif
#if (SPL<=0 && Read_files==1)  // Enforcing
                    if(tiempo_act<d_tCC[timeCC-2]){
                        interpolate_time(d_CC,d_tCC,timeCC,tiempo_act,h1,zb1,q1,npx,npy,npx+npy+i);
                    }
                    else{ // Salida libre
                        q1.y = q0.y;
                    }
#endif

// h1= 4.0*cos(2*PI*tiempo_act/5.0);
// q1.y = h1*q1.y/h0;
                        pos0 *=4;
                        procesarArista_Volumenes(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FD,FU,0,1);
                        
                    }
                    else{
                        h0 = reconsUh[pos0];
                        zb0 = reconsUzb[pos0];
                        q0.x = reconsUhu[pos0].x;
                        q0.y = reconsUhu[pos0].y;
                        
                        h1 = reconsDh[pos1];
                        zb1 = reconsDzb[pos1];
                        q1.x = reconsDhu[pos1].x;
                        q1.y = reconsDhu[pos1].y;
                        
                        
                        pos0 *=4;
                        pos1*=4;
                        procesarArista_Volumenes(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FD,FU,0,1);
                    }
                }
            }
        
    } 
    
    
    
    
}



inline void aristas_volumenes(double *d_dVh, double2 *d_dVhu, double *d_dVH, double *d_Ah, double2 *d_Ahu, double *d_Azb, double *d0_dVh, double2 *d0_dVhu, double *d0_dVH, int npx, int npy, double delta_T, double dx, double dy, double tiempo_act, double epsilon_h, double vmax, double* reconsLh, double* reconsRh, double* reconsDh, double* reconsUh, double* reconsLzb, double* reconsRzb, double* reconsDzb, double* reconsUzb, double2* reconsLhu, double2* reconsRhu, double2* reconsDhu, double2* reconsUhu, double *FL, double *FR, double *FU, double *FD, double *d_CC, double *d_tCC, long unsigned int timeCC)
{
    
    
    
    // Tamano del grid y de bloque en el procesamiento de aristas
//     dim3 threadBlock(NUM_HEBRAS_ANCHO_ARI, NUM_HEBRAS_ALTO_ARI);
    dim3 blockGridEst(iDivUp(npx+1, NUM_HEBRAS_ANCHO_EST), iDivUp(npy+1, NUM_HEBRAS_ALTO_EST));
    dim3 threadBlockEst(NUM_HEBRAS_ANCHO_EST, NUM_HEBRAS_ALTO_EST);
    
    
    volumenes<<<blockGridEst, threadBlockEst>>>(d_dVh, d_dVhu, d_dVH, d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, npx, npy, delta_T, dx, dy, tiempo_act, epsilon_h, vmax, reconsLh, reconsRh, reconsDh, reconsUh, reconsLzb, reconsRzb, reconsDzb, reconsUzb, reconsLhu, reconsRhu, reconsDhu, reconsUhu, FL, FR, FU, FD, d_CC, d_tCC, timeCC);

//     sum_contribuciones<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, npx, npy, delta_T, dx, dy, epsilon_h, vmax, FL, FR, FU, FD);
    
}
#endif 

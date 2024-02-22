

#ifndef _VOLUMEN_KERNEL_H_
#define _VOLUMEN_KERNEL_H_

#include <stdio.h>
#include "Constantes.hxx"
#include "Matriz.cu"



__global__ void obtenerDeltaTVolumenesGPU(double* d_dVh, double2* d_dVhu, double *d_dVH, double* d_deltaTV, int npx, int npy, double dx, double dy, double CFL, double epsilon_h)
{
    int pos,i,j;
    double2 q;
    double2 u;
    double h;
    double cg,lx,ly,dt;
    

    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;
    

    if ((i<npx) && (j<npy)) {
        
        pos=j*npx+i;
        h=d_dVh[pos];
        
        q=d_dVhu[pos];
        calculo_u(h, q, epsilon_h, u);
        
        cg=sqrt(G*h);
        lx=fabs(u.x)+cg;
        ly=fabs(u.y)+cg;
//         dt=CFL*0.5*dx*dy/(dy*lx+dx*ly+EPSILON);
        dt = CFL * 0.5 * min(dx/lx, dy/ly); 
        
        d_deltaTV[pos]=dt;
    }
}

__device__ void interpolate_timeSPL(double *CC, double *tCC, long unsigned int timeCC, double tiempo_act, double &h1, double &zb1, double2 &q1, int npx, int npy, int pos){
    
    
    double dt;
    dt = tCC[1]-tCC[0];
    int t, tp;
    t = int((tiempo_act-tCC[0])/dt);
    t*=(t>0);
    if(t>timeCC-1)
        t=timeCC-1;
    tp = t+1;
    int posCC = t*(2*npx+2*npy)*4 + pos*4;
    int posCCp = tp*(2*npx+2*npy)*4 + pos*4;
    
    double m, b;
    m = (CC[posCCp]-CC[posCC])/(tCC[tp]-tCC[t]);
    b = CC[posCC] - m*tCC[t];
    h1 = m*tiempo_act+b;
    h1*=UNITS;
    
    m = (CC[posCCp+3]-CC[posCC+3])/(tCC[tp]-tCC[t]);
    b = CC[posCC+3] - m*tCC[t];
    zb1 = m*tiempo_act+b;

    h1 += zb1;    

    m = (CC[posCCp+1]-CC[posCC+1])/(tCC[tp]-tCC[t]);
    b = CC[posCC+1] - m*tCC[t];
    q1.x = m*tiempo_act+b;    
    q1.x*= UNITS* h1;

    m = (CC[posCCp+2]-CC[posCC+2])/(tCC[tp]-tCC[t]);
    b = CC[posCC+2] - m*tCC[t];
    q1.y = m*tiempo_act+b;   
    q1.y*= UNITS*h1;

    
}





__global__ void calculo_hmax(double* d_hmax, double * d_dVh, int npx, int npy) {
    int i,j;
    int pos;
    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
    j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;
    if ((i<npx) && (j<npy)) {
       pos=j*npx+i;
       if (d_dVh[pos]>EPSILON) {
         d_hmax[pos]=max(d_hmax[pos],d_dVh[pos]);
       } 
    }

}

__global__ void set_inicializacion(double* d_hmax, int npx, int npy) {
    int i,j;
    int pos;
    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
    j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;
    if ((i<npx) && (j<npy)) {
       pos=j*npx+i;
       d_hmax[pos] = 0.0;
    }
}




__device__ void spongeLayer(double *d_Ah, double2 *d_Ahu, double *d_Azb, int npx, int npy, double *d, int pos, int pos1, int pos2, int posd1, int posd2, double *CC, double *tCC, long unsigned int timeCC, double tiempo_act){
    
    double h1, h2, H1, H2;
    double2 q1, q2;
    

    //     H1=H2=0.0;
    //     h1=h2=1.0;
    //     q1.x=q2.x=0.0;
    //     q1.y=q2.y=0.0;
        
    // Esquina
    if(pos2>0){
        interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, h1, H1, q1, npx, npy, pos1);
        interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, h2, H2, q2, npx, npy, pos2);
        double aux = 0.5*(d[posd1]+d[posd2])*0.5;
        d_Ah[pos] = aux*(h1+h2) + (1.0-d[posd1])*d_Ah[pos]*0.5 + (1.0-d[posd2])*d_Ah[pos] *0.5 ;
        d_Ahu[pos].x = aux*(q1.x+q2.x) + (1.0-d[posd1])*d_Ahu[pos].x*0.5 + (1.0-d[posd2])*d_Ahu[pos].x *0.5 ;
        d_Ahu[pos].y = aux*(q1.y+q2.y) + (1.0-d[posd1])*d_Ahu[pos].y*0.5 + (1.0-d[posd2])*d_Ahu[pos].y *0.5 ;
        d_Azb[pos] = aux*(H1+H2) + (1.0-d[posd1])*d_Azb[pos]*0.5 + (1.0-d[posd2])*d_Azb[pos] *0.5 ;
    }
    // Franja
    else{
        interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, h1, H1, q1, npx, npy, pos1);
        d_Ah[pos]    = d[posd1]*h1 + (1.0-d[posd1])*d_Ah[pos];
        d_Ahu[pos].x = d[posd1]*q1.x + (1.0-d[posd1])*d_Ahu[pos].x;
        d_Ahu[pos].y = d[posd1]*q1.y + (1.0-d[posd1])*d_Ahu[pos].y;
        d_Azb[pos]   = d[posd1]*H1 + (1.0-d[posd1])*d_Azb[pos];
    }


}



__global__ void actualiza_estados(double* d_Ah, double2* d_Ahu, double *d_Azb, double* d0_dVh, double2* d0_dVhu, double *d0_dVH, double* d_dVh, double2* d_dVhu, double *d_dVH, int npx, int npy, int k, double epsilon_h, double invx, double invy, double *FL, double *FR, double *FU, double *FD, double *d, double tiempo_act, double *d_CC, double *d_tCC, long unsigned int timeCC)
{
    int i,j;
    int pos;
    double coef0,coef1;
    double2 q;
    double h;

    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;

    if ((i<npx) && (j<npy)) {
        
        pos = j*npx + i;
        int posF = pos*4;
        
        d_Ah[pos] -= invx*(FL[posF]+FR[posF])  + invy*(FD[posF]+FU[posF]);
        d_Ahu[pos].x -= invx*(FL[posF+1]+FR[posF+1])  + invy*(FD[posF+1]+FU[posF+1]);
        d_Ahu[pos].y -= invx*(FL[posF+2]+FR[posF+2]) + invy*(FD[posF+2]+FU[posF+2]);
        d_Azb[pos] -= invx*(FL[posF+3]+FR[posF+3])  + invy*(FD[posF+3]+FU[posF+3]);
        
        
#if SPL>0        
        if(tiempo_act<d_tCC[timeCC-2]){
            // Franja inferior
            if(j<SPL){
                // Esquina inf izq
                if(i<SPL){
                    spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, i, 2*npx+npy+j, i, j, d_CC, d_tCC, timeCC, tiempo_act);
                }
                // Esquina inf drch 
                else if(i>npx-SPL){
                    spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, i, npx+j, npx-1-i, j, d_CC, d_tCC, timeCC, tiempo_act);
                }
                // Franja central 
                else{
                    spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, i, -1, j, -1, d_CC, d_tCC, timeCC, tiempo_act);
                }
            }
            else{
                // Franja superior
                if(j>npy-SPL){
                    // Esquina superior izq 
                    if(i<SPL){
                        spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, npx+npy+i, 2*npx+npy+j, i, npy-1-j, d_CC, d_tCC, timeCC, tiempo_act);
                    }
                    // Esquina superior drch 
                    else if(i>npx-SPL){
                        spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, npx+npy+i, npx+j, npx-1-i, npy-1-j, d_CC, d_tCC, timeCC, tiempo_act);
                    }
                    // Fraja central
                    else{
                        spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, npx+npy+i, -1, npy-1-j, -1, d_CC, d_tCC, timeCC, tiempo_act);
                    }
                }
                else{
                    // Franja lateral izq
                    if(i<SPL){
                        spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, 2*npx+npy+j, -1, i, -1, d_CC, d_tCC, timeCC, tiempo_act);
                    }
                    // Franja lateral drhc
                    if(i>npx-SPL){
                        spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, npx+j, -1, npx-1-i, -1, d_CC, d_tCC, timeCC, tiempo_act);
                    }    
                }
            }
        }
#endif    
 



#if RK == 1
        
        coef0=0.0*(k==0) + 0.5*((k==1)) ;
        coef1=1.0-coef0;
        
        h   = coef0 * d0_dVh[pos] + coef1 * d_Ah[pos];
        q.x = coef0 * d0_dVhu[pos].x + coef1 * d_Ahu[pos].x;
        q.y = coef0 * d0_dVhu[pos].y + coef1 * d_Ahu[pos].y;
        d_Azb[pos] = coef0 * d0_dVH[pos] + coef1 * d_Azb[pos];
        
        filtro_estado(h,q,epsilon_h);
        
        d_Ah[pos] = h;
        d_Ahu[pos].x = q.x;
        d_Ahu[pos].y = q.y;
        
//          if(j==npy-1) printf("%d   %.15f    %.15f    %.15f   \n",i,d_Ah[pos],d_Ahu[pos].x,d_Ahu[pos].y);
        
        if(k==order-1){
            d0_dVh[pos] = h;
            d0_dVhu[pos].x =q.x;
            d0_dVhu[pos].y =q.y;
            d0_dVH[pos] =d_Azb[pos];
           
            d_dVh[pos] = h;
            d_dVhu[pos].x =q.x;
            d_dVhu[pos].y =q.y;
            d_dVH[pos] =d_Azb[pos];
        }
        else{
            d_dVh[pos] = h;
            d_dVhu[pos].x =q.x;
            d_dVhu[pos].y =q.y;
            d_dVH[pos] =d_Azb[pos];
        }
#else
 
        d_dVh[pos] = d0_dVh[pos] = d_Ah[pos];
        d_dVH[pos] = d0_dVH[pos] = d_Azb[pos];
        d_dVhu[pos].x = d0_dVhu[pos].x = d_Ahu[pos].x;
        d_dVhu[pos].y = d0_dVhu[pos].y = d_Ahu[pos].y;
#endif
        
    }
}





#endif



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
    t = int(tiempo_act/dt);
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


__device__ void sponge_layerBottomAndTop(double* d_Ah, double2* d_Ahu, double *d_Azb, int npx, int npy, double *d, double *CC, double *tCC, long unsigned int timeCC, double tiempo_act, int i, int j) {
    
    int posb, post; 
    posb = j*npx + i;
    post = npx*npy - i -j*npx -1;
    

    double hb, hl, ht, hr;
    double2 qb, qt, qr, ql;
    double Hl, Hr, Ht, Hb;
    

    
    // Bottom       
        // Esquina inf izq y superior izq
        if(i<SPL){
            interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, hb, Hb, qb, npx, npy, i);
            interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, hl, Hl, ql, npx, npy, 2*npx+npy+j);
            
            
            double aux = 0.5*(d[j]+d[i])*0.5;
            
            d_Ah[posb] = aux*(hb+hl) + (1.0-d[i])*d_Ah[posb]*0.5 + (1.0-d[j])*d_Ah[posb] *0.5 ;
            d_Ahu[posb].x = aux*(qb.x+ql.x) + (1.0-d[i])*d_Ahu[posb].x*0.5 + (1.0-d[j])*d_Ahu[posb].x *0.5 ;
            d_Ahu[posb].y = aux*(qb.y+ql.y) + (1.0-d[i])*d_Ahu[posb].y*0.5 + (1.0-d[j])*d_Ahu[posb].y *0.5 ;
            d_Azb[posb] = aux*(Hb+Hl) + (1.0-d[i])*d_Azb[posb]*0.5 + (1.0-d[j])*d_Azb[posb] *0.5 ;
            
            interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, ht, Ht, qt, npx, npy, npx+npy+i);
            interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, hl, Hl, ql, npx, npy, 2*npx+npy*2-j-1);
            
            post =  npx*npy - (npx) - j*npx + i;
            
            d_Ah[post] = aux*(ht+hl) + (1.0-d[i])*d_Ah[post]*0.5 + (1.0-d[j])*d_Ah[post] *0.5 ;
            d_Ahu[post].x = aux*(qt.x+ql.x) + (1.0-d[i])*d_Ahu[post].x*0.5 + (1.0-d[j])*d_Ahu[post].x *0.5 ;
            d_Ahu[post].y = aux*(qt.y+ql.y) + (1.0-d[i])*d_Ahu[post].y*0.5 + (1.0-d[j])*d_Ahu[post].y *0.5 ;
            d_Azb[post] = aux*(Ht+Hl) + (1.0-d[i])*d_Azb[post]*0.5 + (1.0-d[j])*d_Azb[post] *0.5 ;
        }
        
        // Esquina inf drch y superior drch
        else if(i>npx-SPL){
            
            interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, hb, Hb, qb, npx, npy, i);
            interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, hr, Hr, qr, npx, npy, npx+j);

            double aux = 0.5*(d[j]+d[npx-1-i])*0.5;

            
            d_Ah[posb] = aux*(hb+hr) + (1.0-d[npx-1-i])*d_Ah[posb]*0.5 + (1.0-d[j])*d_Ah[posb] *0.5 ;
            d_Ahu[posb].x = aux*(qb.x+qr.x) + (1.0-d[npx-1-i])*d_Ahu[posb].x*0.5 + (1.0-d[j])*d_Ahu[posb].x *0.5 ;
            d_Ahu[posb].y = aux*(qb.y+qr.y) + (1.0-d[npx-1-i])*d_Ahu[posb].y*0.5 + (1.0-d[j])*d_Ahu[posb].y *0.5 ;
            d_Azb[posb] = aux*(Hb+Hr) + (1.0-d[npx-1-i])*d_Azb[posb]*0.5 + (1.0-d[j])*d_Azb[posb] *0.5 ;
            
            interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, ht, Ht, qt, npx, npy, npx+npy+i);
            interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, hr, Hr, qr, npx, npy, npx+npy-1-j);
            
            post = (npy-1)*npx - j*npx + i;
            
            d_Ah[post] = aux*(ht+hr) + (1.0-d[npx-1-i])*d_Ah[post]*0.5 + (1.0-d[j])*d_Ah[post] *0.5 ;
            d_Ahu[post].x = aux*(qt.x+qr.x) + (1.0-d[npx-1-i])*d_Ahu[post].x*0.5 + (1.0-d[j])*d_Ahu[post].x *0.5 ;
            d_Ahu[post].y = aux*(qt.y+qr.y) + (1.0-d[npx-1-i])*d_Ahu[post].y*0.5 + (1.0-d[j])*d_Ahu[post].y *0.5 ;
            d_Azb[post] = aux*(Ht+Hr) + (1.0-d[npx-1-i])*d_Azb[post]*0.5 + (1.0-d[j])*d_Azb[post] *0.5 ;
            
        }
        // Franja central inferior y superior
        else{
            
            interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, hb, Hb, qb, npx, npy, i);
            interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, ht, Ht, qt, npx, npy, npx+npy+i);

            d_Ah[posb]    = d[j]*hb + (1.0-d[j])*d_Ah[posb];
            d_Ahu[posb].x = d[j]*qb.x + (1.0-d[j])*d_Ahu[posb].x;
            d_Ahu[posb].y = d[j]*qb.y + (1.0-d[j])*d_Ahu[posb].y;
            d_Azb[posb]   = d[j]*Hb + (1.0-d[j])*d_Azb[posb];
            
            
            post = (npy-1)*npx - j*npx +i;
            
//             d_Ah[post]    = d[j]*ht + (1.0-d[j])*d_Ah[post];
//             d_Ahu[post].x = d[j]*qt.x + (1.0-d[j])*d_Ahu[post].x;
//             d_Ahu[post].y = d[j]*qt.y + (1.0-d[j])*d_Ahu[post].y;
//             d_Azb[post]   = d[j]*Ht + (1.0-d[j])*d_Azb[post];
            
        }

    
}



__device__ void sponge_layerMid(double* d_Ah, double2* d_Ahu, double *d_Azb, int npx, int npy, double *d, double *CC, double *tCC, long unsigned int timeCC, double tiempo_act, int i, int j) {
    
    int posL, posR; 
    double hr, Hr, hl, Hl;
    double2 qr, ql;
    
    interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, hr, Hr, qr, npx, npy, npx+j);
    interpolate_timeSPL(CC, tCC, timeCC, tiempo_act, hl, Hl, ql, npx, npy, 2*npx+npy+j);
        

//         posL = SPL*npx + j + i*npx;
//         posR = SPL*npx + npx-1 - j + i*npx;
        posL = j*npx + i;
        posR = npx-1 +j*npx - i;
        
        d_Ah[posL]    = d[i]*hl + (1.0-d[i])*d_Ah[posL];
        d_Ahu[posL].x = d[i]*ql.x + (1.0-d[i])*d_Ahu[posL].x;
        d_Ahu[posL].y = d[i]*ql.y + (1.0-d[i])*d_Ahu[posL].y;
        d_Azb[posL]   = d[i]*Hl + (1.0-d[i])*d_Azb[posL];
        
        d_Ah[posR]    = d[i]*hr + (1.0-d[i])*d_Ah[posR];
        d_Ahu[posR].x = d[i]*qr.x + (1.0-d[i])*d_Ahu[posR].x;
        d_Ahu[posR].y = d[i]*qr.y + (1.0-d[i])*d_Ahu[posR].y;
        d_Azb[posR]   = d[i]*Hr + (1.0-d[i])*d_Azb[posR];

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
        if(j<SPL){
            sponge_layerBottomAndTop(d_Ah, d_Ahu, d_Azb, npx, npy, d, d_CC,d_tCC,timeCC,tiempo_act, i,j);
        }
        if(i<SPL && j<npy-SPL && j>SPL){
            sponge_layerMid(d_Ah, d_Ahu, d_Azb, npx, npy, d, d_CC,d_tCC,timeCC, tiempo_act, i,j);    
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

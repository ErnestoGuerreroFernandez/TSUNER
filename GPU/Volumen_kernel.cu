

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
        

        
#if solver==0
        q=d_dVhu[pos];
        calculo_u(h, q, epsilon_h, u);
        
        cg=sqrt(G*h);
        lx=fabs(u.x)+cg;
        ly=fabs(u.y)+cg;
//         dt=CFL*0.5*dx*dy/(dy*lx+dx*ly+EPSILON);
        dt = CFL * 0.5 * min(dx/lx, dy/ly); 
#else
      printf("El calculo de autvalores para el dt con solver>1 no implementado. Fin del programa.");
      exit(-1);
#endif
      
//         printf("i = %d   %.15f  \n",j,dt);
        
        d_deltaTV[pos]=dt;
    }
}



__global__ void sponge_layerBottomAndTop(double* d_Ah, double2* d_Ahu, double *d_Azb, int npx, int npy, double *d, double hl, double Hl, double2 ql, double hr, double Hr, double2 qr, double ht, double Ht, double2 qt, double hb, double Hb, double2 qb, double t) {
    
    int i,j;
    int posb, post; 
    
    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;
    
//     printf("%d \n",i);     
//                     ht  =5.1*cos(2*PI*t/50.0);
//                 ht *=(ht>0.0);
//                 ht += 1.0;
                    ht  =3.0*cos(2*PI*t/5.0);
                ht *=(ht>0.0);
                ht += 1.0;
                
    if ((i<npx) && (j<SPL)) {

        // Bottom
        posb = j*npx + i;
        post = npx*npy - i -j*npx;
//         printf("%.15f \n",d_Ah[posb]);        
        // Esquina inf izq y superior izq
        if(i<SPL){
            
            double aux = 0.5*(d[j]+d[i])*0.5;
            
            d_Ah[posb] = aux*(hb+hl) + (1.0-d[i])*d_Ah[posb]*0.5 + (1.0-d[j])*d_Ah[posb] *0.5 ;
            d_Ahu[posb].x = aux*(qb.x+ql.x) + (1.0-d[i])*d_Ahu[posb].x*0.5 + (1.0-d[j])*d_Ahu[posb].x *0.5 ;
            d_Ahu[posb].y = aux*(qb.y+ql.y) + (1.0-d[i])*d_Ahu[posb].y*0.5 + (1.0-d[j])*d_Ahu[posb].y *0.5 ;
            d_Azb[posb] = aux*(Hb+Hl) + (1.0-d[i])*d_Azb[posb]*0.5 + (1.0-d[j])*d_Azb[posb] *0.5 ;
            
            d_Ah[post] = aux*(ht+hl) + (1.0-d[i])*d_Ah[post]*0.5 + (1.0-d[j])*d_Ah[post] *0.5 ;
            d_Ahu[post].x = aux*(qt.x+ql.x) + (1.0-d[i])*d_Ahu[post].x*0.5 + (1.0-d[j])*d_Ahu[post].x *0.5 ;
            d_Ahu[post].y = aux*(qt.y+ql.y) + (1.0-d[i])*d_Ahu[post].y*0.5 + (1.0-d[j])*d_Ahu[post].y *0.5 ;
            d_Azb[post] = aux*(Ht+Hl) + (1.0-d[i])*d_Azb[post]*0.5 + (1.0-d[j])*d_Azb[post] *0.5 ;
        }
        
        // Esquina inf drch y superior drch
        else if(i>npx-SPL){

            double aux = 0.5*(d[j]+d[npx-1-i])*0.5;
            
            d_Ah[posb] = aux*(hb+hr) + (1.0-d[npx-1-i])*d_Ah[posb]*0.5 + (1.0-d[j])*d_Ah[posb] *0.5 ;
            d_Ahu[posb].x = aux*(qb.x+qr.x) + (1.0-d[npx-1-i])*d_Ahu[posb].x*0.5 + (1.0-d[j])*d_Ahu[posb].x *0.5 ;
            d_Ahu[posb].y = aux*(qb.y+qr.y) + (1.0-d[npx-1-i])*d_Ahu[posb].y*0.5 + (1.0-d[j])*d_Ahu[posb].y *0.5 ;
            d_Azb[posb] = aux*(Hb+Hr) + (1.0-d[npx-1-i])*d_Azb[posb]*0.5 + (1.0-d[j])*d_Azb[posb] *0.5 ;
            
            
            d_Ah[post] = aux*(ht+hr) + (1.0-d[npx-1-i])*d_Ah[post]*0.5 + (1.0-d[j])*d_Ah[post] *0.5 ;
            d_Ahu[post].x = aux*(qt.x+qr.x) + (1.0-d[npx-1-i])*d_Ahu[post].x*0.5 + (1.0-d[j])*d_Ahu[post].x *0.5 ;
            d_Ahu[post].y = aux*(qt.y+qr.y) + (1.0-d[npx-1-i])*d_Ahu[post].y*0.5 + (1.0-d[j])*d_Ahu[post].y *0.5 ;
            d_Azb[post] = aux*(Ht+Hr) + (1.0-d[npx-1-i])*d_Azb[post]*0.5 + (1.0-d[j])*d_Azb[post] *0.5 ;
            
        }
        // Franja central inferior y superior
        else{

            d_Ah[posb]    = d[j]*hb + (1.0-d[j])*d_Ah[posb];
            d_Ahu[posb].x = d[j]*qb.x + (1.0-d[j])*d_Ahu[posb].x;
            d_Ahu[posb].y = d[j]*qb.y + (1.0-d[j])*d_Ahu[posb].y;
            d_Azb[posb]   = d[j]*Hb + (1.0-d[j])*d_Azb[posb];
            
            d_Ah[post]    = d[j]*ht + (1.0-d[j])*d_Ah[post];
            d_Ahu[post].x = d[j]*qt.x + (1.0-d[j])*d_Ahu[post].x;
            d_Ahu[post].y = d[j]*qt.y + (1.0-d[j])*d_Ahu[post].y;
            d_Azb[post]   = d[j]*Ht + (1.0-d[j])*d_Azb[post];
            
        }

    }
}



__global__ void sponge_layerMid(double* d_Ah, double2* d_Ahu, double *d_Azb, int npx, int npy, double *d, double hl, double Hl, double2 ql, double hr, double Hr, double2 qr, double t) {
    
    int i,j;
    int posL, posR; 
    
    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;
    
//                 ht  =5.1*cos(2*PI*t/50.0);
//                 ht *=(ht>0.0);
    
    if ((i<npy-2*SPL) && (j<SPL)) {
        
        // Bottom
        posL = SPL*npx + j + i*npx;
        posR = SPL*npx + npx-1 - j + i*npx;
        
        d_Ah[posL]    = d[j]*hl + (1.0-d[j])*d_Ah[posL];
        d_Ahu[posL].x = d[j]*ql.x + (1.0-d[j])*d_Ahu[posL].x;
        d_Ahu[posL].y = d[j]*ql.y + (1.0-d[j])*d_Ahu[posL].y;
        d_Azb[posL]   = d[j]*Hl + (1.0-d[j])*d_Azb[posL];
        
        d_Ah[posR]    = d[j]*hr + (1.0-d[j])*d_Ah[posR];
        d_Ahu[posR].x = d[j]*qr.x + (1.0-d[j])*d_Ahu[posR].x;
        d_Ahu[posR].y = d[j]*qr.y + (1.0-d[j])*d_Ahu[posR].y;
        d_Azb[posR]   = d[j]*Hr + (1.0-d[j])*d_Azb[posR];

    }
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


__global__ void aplica_fricciones(double* d_Ah, double2* d_Ahu, double* d0_dVh, double2* d0_dVhu,int npx, int npy, double epsilon_h, double dt, double tiempo) {
    
    int i,j;
    int pos;

    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;
    
    
    if ((i<npx) && (j<npy)) {
        
        pos=j*npx+i;
        double h;
        double2 q, u;
        
        h = d0_dVh[pos];
        q.x = d0_dVhu[pos].x;
        q.y = d0_dVhu[pos].y;
        
        calculo_u(h,q,epsilon_h,u);
        double modu = sqrt(u.x*u.x+u.y*u.y);
        
        double c = dt*d_FR*modu/(h+EPSILON);
        
        d_Ahu[pos].x /= (1.0+c);
        d_Ahu[pos].y /= (1.0+c);
        
//         En fricciones se pone el modulo de la velocidad en c y se pone la misma friccion en x e y o bien se calcula cada friccion por separado con u.x e u,y
    }
}


__global__ void actualiza_estados(double* d_Ah, double2* d_Ahu, double *d_Azb, double* d0_dVh, double2* d0_dVhu, double *d0_dVH, double* d_dVh, double2* d_dVhu, double *d_dVH, int npx, int npy, int k, double epsilon_h, double vmax, double dt)
{
    int i,j;
    int pos;
    double coef0,coef1;
    double2 q;
    double h;

    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;

    if ((i<npx) && (j<npy)) {
        pos=(j*npx+i);
        
#if RK == 1
        coef0=0.0*(k==0) + 0.5*((k==1)) ;
        coef1=1.0-coef0;
        
        h   = coef0 * d0_dVh[pos] + coef1 * d_Ah[pos];
        q.x = coef0 * d0_dVhu[pos].x + coef1 * d_Ahu[pos].x;
        q.y = coef0 * d0_dVhu[pos].y + coef1 * d_Ahu[pos].y;
        d_Azb[pos] = coef0 * d0_dVH[pos] + coef1 * d_Azb[pos];
        
        filtro_estado(h,q,epsilon_h,vmax);
        
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


__global__ void borrar(double* d_Ah, double2* d_Ahu, double *d_Azb, double* d0_dVh, double2* d0_dVhu, double *d0_dVH, double* d_dVh, double2* d_dVhu, double *d_dVH, int npx, int npy, int k, double epsilon_h, double vmax, double dt)
{
    int i,j;
    int pos;

    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;

    if ((i<npx) && (j<npy)) {
        pos = j*npx + i;
        
        d_Ah[pos] += 1;
        
        
//                                     int pos0 = (j-1)*npx + i;
//             int pos1 = j*npx + i;
//             if(j==0) pos0=(npy-1)*npx +i;
//             printf("%.15f \n", d0_dVh[pos0]-d0_dVh[pos1]);
        
    }
}



__global__ void borrar2(double* d_Ah, double2* d_Ahu, double *d_Azb, double* d0_dVh, double2* d0_dVhu, double *d0_dVH, double* d_dVh, double2* d_dVhu, double *d_dVH, int npx, int npy, int k, double epsilon_h, double vmax, double dt)
{
    int i,j;
//     int pos;

    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;

    if ((i<npx) && (j<npy)) {
//         pos = j*npx + i;
        
//         d_Ah[pos] += 1;
        
                            int pos0 = (j-1)*npx + i;
            int pos1 = j*npx + i;
            if(j==0) pos0=(npy-1)*npx +i;
            double aux;
            aux = d_Ah[pos0]-d_Ah[pos1];
            if(abs(aux)>1e-14){
            printf("%d    %d   %.15f \n",i, j, aux);    
            }
            
        
    }
}



#endif

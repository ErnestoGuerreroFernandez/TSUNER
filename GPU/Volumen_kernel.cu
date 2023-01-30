

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
        

//         printf("i = %d   %.15f  \n",j,d_dVhu[pos].y);
        
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
        
        d_deltaTV[pos]=dt;
    }
}


__global__ void sponge_layer(double* d_Ah, double2* d_Ahu, double* d0_dVh, double2* d0_dVhu,int npar, double epsilon_h, double dt, double tg1, double tg2, double tg3,  double ca,double H, double Dh, double tiempo) {
    
    int i,j;
    int pos0, pos1; // 0=Izq o abajo. 1=drch o arriba. 
    
    i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;
    
    if ((i<npar) && (j<d_SPL)) {
        
        pos0 = i*npar + j;
        pos1 = i*npar + npar -1 - d_SPL + j;
        

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
        
//         if(i==200) printf("%.15f  \n",d0_dVH[pos]);
        
        if(k==order){
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

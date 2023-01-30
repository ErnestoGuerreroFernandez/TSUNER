#ifndef _RECONSTRUYE_KERNEL_H_
#define _RECONSTRUYE_KERNEL_H_

#include <stdio.h>
#include "Constantes.hxx"
#include "Matriz.cu"

__device__ double avg(double a, double b)
{
    double d,aa,bb;
    aa=fabs(a);
    bb=fabs(b);
    if ((aa+bb) > EPSILON) {
        d=(a*bb+b*aa)/(aa+bb);
    }
    else {
        d=0.0;
    }
    return d;
}

__device__ double minmod(double a, double b, double c)
{
    double d;
    if ((a>0.0) && (b>0.0) && (c>0.0)) {
        d=min(a,min(b,c));
    }
    else {
        if ((a<0.0) && (b<0.0) && (c<0.0)) {
            d=max(a,max(b,c));
        }
        else {
            d=0.0;
        }
    }
    return d;
}




__global__ void reconstruye(double* d_dVh, double2* d_dVhu, double* d_dVH, double* d_hx, double* d_hy, double2* d_hux, double2* d_huy, double* d_zbx,double* d_zby, int npx, int npy, double delta_T, double inv_dx, double inv_dy, double epsilon_h,double2* d_Ahu) {

    int pos,posl,posr;
    int i,j;
    int im,ip,jm,jp;
    
    double dta,Hl,Hi,Hr;
    double hi,hl,hr;
    double2 ui,ul,ur;
    double dl, dr, du;
    double deta;
    
	i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;

    if ((i<npx) && (j<npy)) {
        pos = j*npx + i;
        
#if order == 1
        d_hx[pos]=0.0;
        d_zbx[pos]=0.0;
        d_hux[pos].x=0.0;
        d_hux[pos].y=0.0;
        
        d_hy[pos]=0.0;
        d_zby[pos]=0.0;
        d_huy[pos].x=0.0;
        d_huy[pos].y=0.0;
#else   
 
        
        // **********************
        // Direccion x
        // **********************

        hi=d_dVh[pos];
        Hi=d_dVH[pos];
        
        im=i-1*(i>0);
        ip=i+1*(i<npx-1);
        posl=j*npx + im;
        posr=j*npx + ip;
        hl=d_dVh[posl];
        Hl=d_dVH[posl];
        hr=d_dVh[posr];
        Hr=d_dVH[posr];
        dta=delta_T*inv_dx;

        
        // Reconstruccion de h
        dl = hi - hl;
        dr = hr - hi;
        d_hx[pos] = avg(dl,dr);
        if(0.5*fabs(d_hx[pos])>hi){
            d_hx[pos] = 2.0*hi*SGN(d_hx[pos]);
        }
        
        // Reconstruccion de eta y zb
        double zm;
        zm = min(Hl,Hi);
        dl = max(hi-Hi+zm,0.0) - max(hl-Hl+zm,0.0);
        zm = min(Hr,Hi);
        dr = max(hr-Hr+zm,0.0) - max(hi-Hi+zm,0.0);
        deta = avg(dl,dr);
        d_zbx[pos] = d_hx[pos] - deta;
        
        // Reconstruccion de u
        calculo_u(hi,d_dVhu[pos],epsilon_h,ui);
        calculo_u(hl,d_dVhu[posl],epsilon_h,ul);
        calculo_u(hr,d_dVhu[posr],epsilon_h,ur);
        dl = ui.x - ul.x;
        dr = ur.x - ui.x;
        du = avg(dl,dr);
        d_hux[pos].x = hi*du + ui.x * d_hx[pos];
        dl = ui.y - ul.y;
        dr = ur.y - ui.y;
        du = avg(dl,dr);
        d_hux[pos].y = hi*du + ui.y * d_hx[pos];
        
        // Integral
        d_Ahu[pos].x -= dta*G*hi*deta;

        
        // **********************
        // Direccion y
        // **********************
        

        jm=j-1*(j>0);
        jp=j+1*(j<npy-1);
        posl=jm*npx +i;
        posr=jp*npx +i;
        hl=d_dVh[posl];
        Hl=d_dVH[posl];
        hr=d_dVh[posr];
        Hr=d_dVH[posr];
        dta=delta_T*inv_dy;
        
        // Reconstruccion de h
        dl = hi - hl;
        dr = hr - hi;
        d_hy[pos] = avg(dl,dr);
        if(0.5*fabs(d_hy[pos])>hi){
            d_hy[pos] = 2.0*hi*SGN(d_hy[pos]);
        }
        
        
        // Reconstruccion de eta y zb
        zm = min(Hl,Hi);
        dl = max(hi-Hi+zm,0.0) - max(hl-Hl+zm,0.0);
        zm = min(Hr,Hi);
        dr = max(hr-Hr+zm,0.0) - max(hi-Hi+zm,0.0);
        deta = avg(dl,dr);
        d_zby[pos] = d_hy[pos] - deta;
        
        
        // Reconstruccion de u
        calculo_u(hl,d_dVhu[posl],epsilon_h,ul);
        calculo_u(hr,d_dVhu[posr],epsilon_h,ur);
        dl = ui.x - ul.x;
        dr = ur.x - ui.x;
        du = avg(dl,dr);
        d_huy[pos].x = hi*du + ui.x * d_hy[pos];
        dl = ui.y - ul.y;
        dr = ur.y - ui.y;
        du = avg(dl,dr);
        d_huy[pos].y = hi*du + ui.y * d_hy[pos];
        
        // Integral
        d_Ahu[pos].y -= dta*G*hi*deta;

        
#endif
        
    }

}

#endif

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

__device__ void aplica_fricciones(double2 &d_Ahu, double d0_dVh, double2 &d0_dVhu, double epsilon_h, double dt) {
    
    double h;
    double2 q, u;
    
    h = d0_dVh;
    q.x = d0_dVhu.x;
    q.y = d0_dVhu.y;
    
    calculo_u(h,q,epsilon_h,u);
    double modu = sqrt(u.x*u.x+u.y*u.y);
    
    double c = dt*d_FR*modu/(h+EPSILON);
    
    d_Ahu.x /= (1.0+c);
    d_Ahu.y /= (1.0+c);
            
}



__global__ void reconstruye(double* d_dVh, double2* d_dVhu, double* d_dVH, double* reconsLh, double* reconsRh, double* reconsDh, double* reconsUh, double* reconsLzb, double* reconsRzb, double* reconsDzb, double* reconsUzb, double2* reconsLhu, double2* reconsRhu, double2* reconsDhu, double2* reconsUhu,int npx, int npy, double delta_T, double inv_dx, double inv_dy, double epsilon_h,double2* d_Ahu, double *d0_dVh, double2 *d0_dVhu) {


    int i,j;
    int pos;
	i = blockIdx.x*NUM_HEBRAS_ANCHO_EST + threadIdx.x;
	j = blockIdx.y*NUM_HEBRAS_ALTO_EST + threadIdx.y;

    if ((i<npx) && (j<npy)) {
        pos = j*npx + i;
        
#if order == 1

        reconsLh[pos] = reconsRh[pos] = reconsUh[pos] = reconsDh[pos] = d_dVh[pos];
        reconsLzb[pos] = reconsRzb[pos] = reconsUzb[pos] = reconsDzb[pos] = d_dVH[pos];
        reconsLhu[pos].x = reconsRhu[pos].x = reconsUhu[pos].x = reconsDhu[pos].x = d_dVhu[pos].x;
        reconsLhu[pos].y = reconsRhu[pos].y = reconsUhu[pos].y = reconsDhu[pos].y = d_dVhu[pos].y;
        
#else   
 
        int posl,posr;
        int im,ip,jm,jp;
        double dta,Hl,Hi,Hr;
        double hi,hl,hr;
        double2 ui,ul,ur;
        double dl, dr, du;
        double deta;
        
        // **********************
        // Direccion x
        // **********************

        hi=d_dVh[pos];
        Hi=d_dVH[pos];
        
        im=i-1*(i>0);
        ip=i+1*(i<npx-1);
#if WBC==2
        if(i==0) im = npx-1;
        if(i==npx-1) ip = 0;
#endif
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
        double dh = avg(dl,dr);
        if(0.5*fabs(dh)>hi){
            dh = 2.0*hi*SGN(dh);
        }
        reconsLh[pos] = hi - 0.5*dh;
        reconsRh[pos] = hi + 0.5*dh;
        
        // Reconstruccion de eta y zb
        double zm;
        zm = min(Hl,Hi);
        dl = max(hi-Hi+zm,0.0) - max(hl-Hl+zm,0.0);
        zm = min(Hr,Hi);
        dr = max(hr-Hr+zm,0.0) - max(hi-Hi+zm,0.0);
        deta = avg(dl,dr);
        double dz = dh - deta;
        
        reconsLzb[pos] = Hi - 0.5*dz;
        reconsRzb[pos] = Hi + 0.5*dz;
        
        // Reconstruccion de u
        calculo_u(hi,d_dVhu[pos],epsilon_h,ui);
        calculo_u(hl,d_dVhu[posl],epsilon_h,ul);
        calculo_u(hr,d_dVhu[posr],epsilon_h,ur);
        dl = ui.x - ul.x;
        dr = ur.x - ui.x;
        du = avg(dl,dr);
        double dq = hi*du + ui.x * dh;
        
        reconsLhu[pos].x = d_dVhu[pos].x - 0.5*dq;
        reconsRhu[pos].x = d_dVhu[pos].x + 0.5*dq;
        
        dl = ui.y - ul.y;
        dr = ur.y - ui.y;
        du = avg(dl,dr);
        dq = hi*du + ui.y * dh;
        
        reconsLhu[pos].y = d_dVhu[pos].y - 0.5*dq;
        reconsRhu[pos].y = d_dVhu[pos].y + 0.5*dq;
        
        // Integral
        d_Ahu[pos].x -= dta*G*hi*deta;
        
#if RK==0
        double Tpx = G*hi*deta;
#endif

        
        // **********************
        // Direccion y
        // **********************
        

        jm=j-1*(j>0);
        jp=j+1*(j<npy-1);
#if WBC==2
        if(j==0) jm = npy-1;
        if(j==npy-1) jp = 0;
#endif
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
        dh = avg(dl,dr);
        if(0.5*fabs(dh)>hi){
            dh = 2.0*hi*SGN(dh);
        }
        reconsDh[pos] = hi - 0.5*dh;
        reconsUh[pos] = hi + 0.5*dh;
        
        
        // Reconstruccion de eta y zb
        zm = min(Hl,Hi);
        dl = max(hi-Hi+zm,0.0) - max(hl-Hl+zm,0.0);
        zm = min(Hr,Hi);
        dr = max(hr-Hr+zm,0.0) - max(hi-Hi+zm,0.0);
        deta = avg(dl,dr);
        dz = dh - deta;
        
        reconsDzb[pos] = Hi - 0.5*dz;
        reconsUzb[pos] = Hi + 0.5*dz;
        

        // Reconstruccion de u
        calculo_u(hl,d_dVhu[posl],epsilon_h,ul);
        calculo_u(hr,d_dVhu[posr],epsilon_h,ur);
        dl = ui.x - ul.x;
        dr = ur.x - ui.x;
        du = avg(dl,dr);
        dq = hi*du + ui.x * dh;
        
        reconsDhu[pos].x  = d_dVhu[pos].x - 0.5*dq;
        reconsUhu[pos].x  = d_dVhu[pos].x + 0.5*dq;
        
        dl = ui.y - ul.y;
        dr = ur.y - ui.y;
        du = avg(dl,dr);
        dq = hi*du + ui.y * dh;
        
        reconsDhu[pos].y  = d_dVhu[pos].y - 0.5*dq;
        reconsUhu[pos].y  = d_dVhu[pos].y + 0.5*dq;
        
        // Integral
        d_Ahu[pos].y -= dta*G*hi*deta;
        
#if RK==0
        double Tpy = G*hi*deta;
#endif
        
                
#if RK==0            
        double2 uRR, uL, uU, uD;
        calculo_u(reconsRh[pos],reconsRhu[pos],epsilon_h,uRR);
        calculo_u(reconsLh[pos],reconsLhu[pos],epsilon_h,uL);
        calculo_u(reconsUh[pos],reconsUhu[pos],epsilon_h,uU);
        calculo_u(reconsDh[pos],reconsDhu[pos],epsilon_h,uD);
    
        double qhx, qhy;
        qhx = (reconsRhu[pos].x-reconsLhu[pos].x) *inv_dx;
        qhy = (reconsUhu[pos].y-reconsDhu[pos].y) *inv_dy;
        
        
        // MUSCL-Hancock
        reconsLh[pos] -= 0.5 * delta_T * (qhx + qhy);
        reconsRh[pos] -= 0.5 * delta_T * (qhx + qhy);
        reconsDh[pos] -= 0.5 * delta_T * (qhx + qhy);
        reconsUh[pos] -= 0.5 * delta_T * (qhx + qhy);
        

        qhx = (reconsRhu[pos].x*uRR.x-reconsLhu[pos].x*uL.x)*inv_dx;
        qhy = (reconsUhu[pos].x*uU.y-reconsDhu[pos].x*uD.y)*inv_dy;
        reconsLhu[pos].x -= 0.5 * delta_T * (qhx + qhy);
        reconsRhu[pos].x -= 0.5 * delta_T * (qhx + qhy);
        reconsDhu[pos].x -= 0.5 * delta_T * (qhx + qhy);
        reconsUhu[pos].x -= 0.5 * delta_T * (qhx + qhy);
        

        qhx = (reconsRhu[pos].x*uRR.y - reconsLhu[pos].x*uL.y) *inv_dx;
        qhy = (reconsUhu[pos].y*uU.y - reconsDhu[pos].y*uD.y) *inv_dy;
        reconsLhu[pos].y -= 0.5 * delta_T * (qhx + qhy);
        reconsRhu[pos].y -= 0.5 * delta_T * (qhx + qhy);
        reconsDhu[pos].y -= 0.5 * delta_T * (qhx + qhy);
        reconsUhu[pos].y -= 0.5 * delta_T * (qhx + qhy);
        
        
        qhx = (bedload(reconsRh[pos],uRR)-bedload(reconsLh[pos],uL))*inv_dx;
        qhy = (bedload(reconsUh[pos],uU)-bedload(reconsDh[pos],uD))*inv_dy;
        reconsLzb[pos] -= 0.5 * delta_T * (qhx + qhy);
        reconsRzb[pos] -= 0.5 * delta_T * (qhx + qhy);
        reconsDzb[pos] -= 0.5 * delta_T * (qhx + qhy);
        reconsUzb[pos] -= 0.5 * delta_T * (qhx + qhy);
        
        
        
        reconsLhu[pos].x -= 0.5 * delta_T*inv_dx * Tpx;
        reconsUhu[pos].x -= 0.5 * delta_T*inv_dx * Tpx;
        reconsDhu[pos].x -= 0.5 * delta_T*inv_dx * Tpx;
        reconsRhu[pos].x -= 0.5 * delta_T*inv_dx * Tpx;
        
        reconsLhu[pos].y -= 0.5 * delta_T*inv_dy * Tpy;
        reconsUhu[pos].y -= 0.5 * delta_T*inv_dy * Tpy;
        reconsDhu[pos].y -= 0.5 * delta_T*inv_dy * Tpy;
        reconsRhu[pos].y -= 0.5 * delta_T*inv_dy * Tpy;

#endif   
        
        
        // Filtros
        filtro_estado(reconsLh[pos],reconsLhu[pos],epsilon_h);
        filtro_estado(reconsRh[pos],reconsRhu[pos],epsilon_h);
        filtro_estado(reconsUh[pos],reconsUhu[pos],epsilon_h);
        filtro_estado(reconsDh[pos],reconsDhu[pos],epsilon_h);
        
        
#endif
        // Fricciones
        aplica_fricciones(d_Ahu[pos], d0_dVh[pos], d0_dVhu[pos], epsilon_h, delta_T);        
    }

}

#endif

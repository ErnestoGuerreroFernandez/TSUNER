#ifndef _ARISTA_KERNEL_H_
#define _ARISTA_KERNEL_H_

#include <stdio.h>
#include "Matriz.cu"
#include "Constantes.hxx"







__device__ void procesarArista(double hhl, double zb0, double2 qql,  double hhr, double zb1,double2 qqr,int normal_x, int normal_y, double dl, double dt, double* d_Ah, double2* d_Ahu, double *d_Azb, int pos_vol0, int pos_vol1, double epsilon_h, int frontera)
{
    double2 q0,q1;
    double2 ql, qr;
    double2 u0,u1;
    double h0, h1;
    double dtaux=dt/dl;
    

    
    // x = normal, y = tangencial
    ql.x = qql.x*normal_x+qql.y*normal_y;
    ql.y = qql.y*normal_x-qql.x*normal_y;
    
    qr.x = qqr.x*normal_x+qqr.y*normal_y;
    qr.y = qqr.y*normal_x-qqr.x*normal_y;
    
    calculo_u(hhl,ql,epsilon_h,u0);
    calculo_u(hhr,qr,epsilon_h,u1);
    
    // Reconstruccion hidrostatica
    double zmH = min(zb0,zb1);
    h0 = max(hhl-zb0+zmH,0.0);
    h1 = max(hhr-zb1+zmH,0.0);
    
    q0.x = u0.x * h0;
    q0.y = u0.y * h0;
    
    q1.x = u1.x * h1;
    q1.y = u1.y * h1;
    
    double deta = h1-h0;

    double F0[4], F1[4];
    F0[0] = F0[1] = F0[2] = F0[3] = 0.0;
    F1[0] = F1[1] = F1[2] = F1[3] = 0.0;
    flujoC(h0,u0,F0,1, zb1-zb0, dl);  // Apunte: Puesto que el sistema de referencia es con la H medida por encima de la batimetria, el sistema que resuelves es dt(H) - dx(qb) (ver articulo de Castro-Fernandez-Ferreiro). Por tanto, se pone un -1 en el flujo del sedimento.
    flujoC(h1,u1,F1,1, zb1-zb0, dl);

    double hm, zm; // Medias de Roe
    double2 um;
    roe_average(h0,u0,zb0,h1,u1,zb1,hm,um,zm);
    
    double autovalores[3], SL, SR;
    
    calc_autoval_aprox0(h0,u0,autovalores);
    SL = autovalores[0];
    
    calc_autoval_aprox0(h1,u1,autovalores);
    SR = autovalores[2];
    
    calc_autoval_aprox0(hm,um,autovalores);
    SL=min(SL,autovalores[0]);
    SR=max(SR,autovalores[2]);
    
    double a0, a1;
    if (SR-SL>EPSILON) {
        a0=(SR*fabs(SL)-SL*fabs(SR))/(SR-SL);
        a1=(fabs(SR)-fabs(SL))/(SR-SL);
    }
    else {
        a0=max(fabs(SL),fabs(SR));
        a1=0.;
    }        
    
    
    double df[4], vis[4], DU[4];
    double s0,s1,s2;
    
    df[0] = F1[0]-F0[0];
    df[1] = F1[1]-F0[1];
    df[2] = F1[2]-F0[2];
    df[3] = F1[3]-F0[3];
    
    df[1] += G*hm*deta;
    
    DU[0] = h1-h0;
    DU[1] = q1.x-q0.x;
    DU[2] = q1.y-q0.y;
    DU[3] = zb1-zb0;
    
    s0 = SGN(zb1-zb0);
    s1 = SGN(df[3]);
    s2 = SGN(a1);
    
    if (fabs(s0)>EPSILON){
        s1=s1/s0;
    }
    else{
        s1=0.0;
    }
    
    vis[0] = a0*DU[0]+a1*df[0];
    vis[1] = a0*DU[1]+a1*df[1];
    vis[2] = a0*DU[2]+a1*df[2];
    vis[3] = a0*DU[3]+a1*df[3];
    
    if(s2==s1){
        vis[3] = a1*df[3];
    }
    else{
        vis[3] = -a1*df[3];
    }
    

    double Dm[4], Dp[4];
    Dm[0] = Dm[1] = Dm[2] = Dm[3] = 0.0;
    Dp[0] = Dp[1] = Dp[2] = Dp[3] = 0.0;
    for (int i=0; i<4; i++){
        Dm[i] = 0.5*(df[i]-vis[i]) + F0[i];        
        Dp[i] = 0.5*(df[i]+vis[i]) - F1[i];
    }
 
    
    d_Ah[pos_vol0]-=dtaux*Dm[0];

    d_Ahu[pos_vol0].x-=dtaux*(normal_x*Dm[1]-normal_y*Dm[2]);
    d_Ahu[pos_vol0].y-=dtaux*(normal_y*Dm[1]+normal_x*Dm[2]);
    
    
    d_Azb[pos_vol0] -=dtaux*Dm[3];
    
    #if WBC==2
    d_Ah[pos_vol1]-=dtaux*Dp[0];    
    d_Ahu[pos_vol1].x-=dtaux*(normal_x*Dp[1]-normal_y*Dp[2]);
    d_Ahu[pos_vol1].y-=dtaux*(normal_y*Dp[1]+normal_x*Dp[2]);    
    d_Azb[pos_vol1] -= dtaux*Dp[3];    
    
    #else
    
    if (pos_vol1 != -1) {
        
        d_Ah[pos_vol1]-=dtaux*Dp[0];
        
        d_Ahu[pos_vol1].x-=dtaux*(normal_x*Dp[1]-normal_y*Dp[2]);
        d_Ahu[pos_vol1].y-=dtaux*(normal_y*Dp[1]+normal_x*Dp[2]);
        
        d_Azb[pos_vol1] -= dtaux*Dp[3];
    }
    #endif
    

}


__global__ void procesarAristasVerGPU(double* d_dVH, double* d_dVh, double2* d_dVhu, double* d_hx, double2* d_hux,double* d_zbx, int npx, int npy, double delta_T, double* d_Ah, double2* d_Ahu, double *d_Azb, double dx, double tiempo_act,int tipo, double epsilon_h, double vmax, double t)
{
    int i, j;
    double2 ql,qr;
    double hl,hr;
    double zbl,zbr;
    int pos0,pos1;
    int frontera=0;
    pos1=-1;
    

    if(tipo==1){
        i = 2*(blockIdx.x*NUM_HEBRAS_ANCHO_ARI + threadIdx.x);
        j = blockIdx.y*NUM_HEBRAS_ALTO_ARI + threadIdx.y;
    }
    else{
        i = 2*(blockIdx.x*NUM_HEBRAS_ANCHO_ARI + threadIdx.x) + 1;
        j = blockIdx.y*NUM_HEBRAS_ALTO_ARI + threadIdx.y;
    }

    
    if ((i <= npx) && (j < npy)) {
        if (i == 0) {
            #if WBC==2

            #else

            frontera=4;
            // Frontera izquierda
            // normal_x = -1.0
            // El volumen 0 de la arista esta situado a la derecha de la arista
            pos0=npx*j + i;
//             evalua_rh(d_dVh[pos0],d_hx[pos0],-0.5,hl);
            evalua_rzb(d_dVH[pos0],d_zbx[pos0],-0.5,zbl);
            evalua_rq(d_dVhu[pos0],d_hux[pos0],-0.5,ql);
            
            hl = d_dVh[pos0];
            
            filtro_estado(hl, ql , epsilon_h);
            
            hr=hl;
            zbr=zbl;
            qr.y=ql.y;
            #if WBC==1
            qr.x=-ql.x;
            #else
            qr.x=ql.x;
            #endif
            //             hr = 3.0*cos(2*PI*t/5.0);
            //             qr.x = hr * qr.x/(hr+EPSILON);
            //             
            //             qr.y = 0.0;
            
            procesarArista(hl, zbl, ql, hr, zbr, qr,-1, 0, dx, delta_T,d_Ah, d_Ahu, d_Azb, pos0,pos1, epsilon_h,frontera);
            #endif// Del WBC==2
        }
        else {
            // normal_x = 1.0
            // El volumen 0 de la arista esta situado a la izquierda de la arista
            pos0 = npx*j + i-1;
//             evalua_rh(d_dVh[pos0],d_hx[pos0],0.5,hl);
            evalua_rzb(d_dVH[pos0],d_zbx[pos0],0.5,zbl);
            evalua_rq(d_dVhu[pos0],d_hux[pos0],0.5,ql);
            
            hl = d_dVh[pos0];
            
            filtro_estado(hl, ql , epsilon_h);
            
            if (i == npx) {
                frontera=2;
                hr=hl;
                zbr=zbl;
                #if WBC==1
                qr.x=-ql.x;
                #else
                qr.x=ql.x;
                #endif
                qr.y=ql.y;
                
                #if WBC==2
                pos1 = npx*j + 0;
//                 evalua_rh(d_dVh[pos1],d_hx[pos1],-0.5,hr);
                evalua_rzb(d_dVH[pos1],d_zbx[pos1],-0.5,zbr);
                evalua_rq(d_dVhu[pos1],d_hux[pos1],-0.5,qr);
                
                hr = d_dVh[pos1];
                
                filtro_estado(hr, qr , epsilon_h);  
                
                #endif
                    
                procesarArista(hl, zbl, ql, hr, zbr, qr, 1, 0, dx, delta_T,d_Ah, d_Ahu, d_Azb, pos0,pos1, epsilon_h,frontera);
            }
            else {
                // Arista interna
                // El volumen 1 de la arista est� situado a la derecha de la arista
                pos1 = npx*j + i;
//                 evalua_rh(d_dVh[pos1],d_hx[pos1],-0.5,hr);
                evalua_rzb(d_dVH[pos1],d_zbx[pos1],-0.5,zbr);
                evalua_rq(d_dVhu[pos1],d_hux[pos1],-0.5,qr);
                
                hr = d_dVh[pos1];
                
                filtro_estado(hr, qr , epsilon_h);
                procesarArista(hl, zbl, ql, hr, zbr, qr, 1, 0, dx, delta_T,d_Ah, d_Ahu, d_Azb, pos0,pos1, epsilon_h,frontera);
                
            }
        }
    }
}

__global__ void procesarAristasHorGPU(double* d_dVH, double* d_dVh, double2* d_dVhu, double* d_hy, double2* d_huy, double* d_zby,int npx, int npy, double delta_T, double* d_Ah, double2* d_Ahu, double *d_Azb, double dy, double tiempo_act, int tipo,double epsilon_h,double vmax, double t)
{
    int i, j;
    
    double2 ql,qr;
    double hl,hr;
    double zbl,zbr;
    
    int pos0,pos1;
    int frontera=0;
    pos1 = -1;
    
    // Multiplicamos j por 2 porque se procesan aristas alternas
//     i = blockIdx.x*NUM_HEBRAS_ANCHO_ARI + threadIdx.x;
//     j = 2*(blockIdx.y*NUM_HEBRAS_ALTO_ARI + threadIdx.y);
//     if (tipo == 2) j++;
    if(tipo==1){
        i = blockIdx.x*NUM_HEBRAS_ANCHO_ARI + threadIdx.x;
        j = 2*(blockIdx.y*NUM_HEBRAS_ALTO_ARI + threadIdx.y);        
    }
    else{
        i = blockIdx.x*NUM_HEBRAS_ANCHO_ARI + threadIdx.x;
        j = 2*(blockIdx.y*NUM_HEBRAS_ALTO_ARI + threadIdx.y) + 1;
    }
    
    // Comprobamos si la hebra (arista) est� dentro de los l�mites de la malla
    if ((i < npx) && (j <= npy)) {
        // Obtenemos los datos de los vol�menes 0 y 1
        if (j == 0) {
            // Frontera inferior
            // normal_y = -1.0
            // El volumen 0 de la arista esta situado debajo de la arista
            #if WBC==2
            //                 continue;
            #else
            frontera=1;
            pos0 = npx*j + i;
            
//             evalua_rh(d_dVh[pos0],d_hy[pos0],-0.5,hl);
            evalua_rzb(d_dVH[pos0],d_zby[pos0],-0.5,zbl);
            evalua_rq(d_dVhu[pos0],d_huy[pos0],-0.5,ql);
            
            hl = d_dVh[pos0];
            
            filtro_estado(hl, ql , epsilon_h);
            
            hr=hl;
            zbr=zbl;
            qr.x=ql.x;
            #if WBC==1
            qr.y=-ql.y;
            #else
            qr.y=ql.y;
            #endif
            //             hr = 3.0*cos(2*PI*t/5.0);
            //             qr.y = hr * qr.x/(hr+EPSILON);
            //             
            //             qr.x = 0.0;
            pos1=-1;
            procesarArista(hl, zbl, ql, hr, zbr, qr, 0, -1, dy, delta_T,d_Ah, d_Ahu, d_Azb, pos0,pos1, epsilon_h,frontera);
            #endif // Del WBC==2
            
        }
        else {
            // normal_y = 1.0
            // El volumen 0 de la arista esta situado arriba de la arista
            pos0 = npx*(j-1) + i;
//             evalua_rh(d_dVh[pos0],d_hy[pos0],0.5,hl);
            evalua_rzb(d_dVH[pos0],d_zby[pos0],0.5,zbl);
            evalua_rq(d_dVhu[pos0],d_huy[pos0],0.5,ql);
            
            hl = d_dVh[pos0];
            
            filtro_estado(hl, ql , epsilon_h);
            
            if (j == npy) {
                frontera=3;
                // Frontera superior. El volumen 1 es fantasma
                hr=hl;
                zbr=zbl;
                qr.x=ql.x;
                #if WBC==1
                qr.y=-ql.y;
                #else
                qr.y=ql.y;
                #endif
                
                #if WBC==2
                pos1 = i;
//                 evalua_rh(d_dVh[pos1],d_hy[pos1],-0.5,hr);
                evalua_rzb(d_dVH[pos1],d_zby[pos1],-0.5,zbr);
                evalua_rq(d_dVhu[pos1],d_huy[pos1],-0.5,qr);
                
                hr = d_dVh[pos1];
                
                filtro_estado(hr, qr , epsilon_h);
                #endif
                
                pos1=-1;
                procesarArista(hl, zbl, ql, hr, zbr, qr,0, 1, dy, delta_T,d_Ah, d_Ahu, d_Azb, pos0,pos1, epsilon_h,frontera);
            }
            else {
                // Arista interna
                // El volumen 1 de la arista esta situado debajo de la arista
                pos1 = npx*j + i;
//                 evalua_rh(d_dVh[pos1],d_hy[pos1],-0.5,hr);
                evalua_rzb(d_dVH[pos1],d_zby[pos1],-0.5,zbr);
                evalua_rq(d_dVhu[pos1],d_huy[pos1],-0.5,qr);
                
                hr = d_dVh[pos1];
                
                filtro_estado(hr, qr , epsilon_h);
                
//                 if(abs(hr-hl)>1e-11) printf("%d   %d   %.15f  \n", i, j, hr-hl);
                procesarArista(hl, zbl, ql, hr, zbr, qr,0, 1, dy, delta_T,d_Ah, d_Ahu, d_Azb, pos0,pos1, epsilon_h,frontera);
            }
        }
    }
}



__device__ void procesarArista_Volumenes(double hhl, double zb0, double2 qql, double hhr, double zb1, double2 qqr, int pos0, int pos1, double epsilon_h, double vmax, double *FL, double *FR, int normal_x, int normal_y, int fixed_sediment, double dl){
    
    double2 q0,q1;
    double2 ql, qr;
    double2 u0,u1;
    double h0, h1;
    

    
    // x = normal, y = tangencial
    ql.x = qql.x*normal_x+qql.y*normal_y;
    ql.y = qql.y*normal_x-qql.x*normal_y;
    
    qr.x = qqr.x*normal_x+qqr.y*normal_y;
    qr.y = qqr.y*normal_x-qqr.x*normal_y;
    
    calculo_u(hhl,ql,epsilon_h,u0);
    calculo_u(hhr,qr,epsilon_h,u1);
    
    // Reconstruccion hidrostatica
    double zmH = min(zb0,zb1);
    h0 = max(hhl-zb0+zmH,0.0);
    h1 = max(hhr-zb1+zmH,0.0);
    
    q0.x = u0.x * h0;
    q0.y = u0.y * h0;
    
    q1.x = u1.x * h1;
    q1.y = u1.y * h1;
    
    double deta = h1-h0;

    double F0[4], F1[4];
    F0[0] = F0[1] = F0[2] = F0[3] = 0.0;
    F1[0] = F1[1] = F1[2] = F1[3] = 0.0;
    flujoC(h0,u0,F0,fixed_sediment, zb1-zb0, dl);  // Apunte 1: Puesto que el sistema de referencia es con la H medida por encima de la batimetria, el sistema que resuelves es dt(H) - dx(qb) (ver articulo de Castro-Fernandez-Ferreiro). Por tanto, se pone un -1 en el flujo del sedimento.
    flujoC(h1,u1,F1,fixed_sediment, zb1-zb0, dl);

    double hm, zm; // Medias de Roe
    double2 um;
    roe_average(h0,u0,zb0,h1,u1,zb1,hm,um,zm);
    
    double autovalores[3], SL, SR;
    
    calc_autoval_aprox0(h0,u0,autovalores);
    SL = autovalores[0];
    
    calc_autoval_aprox0(h1,u1,autovalores);
    SR = autovalores[2];
    
    calc_autoval_aprox0(hm,um,autovalores);
    SL=min(SL,autovalores[0]);
    SR=max(SR,autovalores[2]);
    
    double a0, a1;
    if (SR-SL>EPSILON) {
        a0=(SR*fabs(SL)-SL*fabs(SR))/(SR-SL);
        a1=(fabs(SR)-fabs(SL))/(SR-SL);
    }
    else {
        a0=max(fabs(SL),fabs(SR));
        a1=0.;
    }        
    
    
    double df[4], vis[4], DU[4];
    double s0,s1,s2;
    
    df[0] = F1[0]-F0[0]  ;
    df[1] = F1[1]-F0[1];
    df[2] = F1[2]-F0[2];
    df[3] = F1[3]-F0[3];
    
    df[1] += G*hm*deta  ;
    
    DU[0] = h1-h0;
    DU[1] = q1.x-q0.x;
    DU[2] = q1.y-q0.y;
    DU[3] = zb1-zb0;
    
    s0 = SGN(zb1-zb0);
    s1 = SGN(df[3]);
    s2 = SGN(a1);
    
    if (fabs(s0)>EPSILON){
        s1=s1/s0;
    }
    else{
        s1=0.0;
    }
    
    vis[0] = a0*DU[0]+a1*df[0];
    vis[1] = a0*DU[1]+a1*df[1];
    vis[2] = a0*DU[2]+a1*df[2];
    vis[3] = a0*DU[3]+a1*df[3];
    
    if(s2==s1){
        vis[3] = a1*df[3];
    }
    else{
        vis[3] = -a1*df[3];
    }
    

    double Dm[4], Dp[4];
    Dm[0] = Dm[1] = Dm[2] = Dm[3] = 0.0;
    Dp[0] = Dp[1] = Dp[2] = Dp[3] = 0.0;
    for (int i=0; i<4; i++){
        Dm[i] = 0.5*(df[i]-vis[i]) + F0[i]  ;
        Dp[i] = 0.5*(df[i]+vis[i]) - F1[i]  ;
    }

 

    
#if WBC==2
    FR[pos0]   = Dm[0];
    FR[pos0+1] = (normal_x*Dm[1]-normal_y*Dm[2]);
    FR[pos0+2] = (normal_y*Dm[1]+normal_x*Dm[2]);
    FR[pos0+3] = Dm[3];
    
    FL[pos1]   = Dp[0];
    FL[pos1+1] = (normal_x*Dp[1]-normal_y*Dp[2]);
    FL[pos1+2] = (normal_y*Dp[1]+normal_x*Dp[2]);
    FL[pos1+3] = Dp[3];
#else
    if (pos1>0){
        FR[pos0]   = Dm[0];
        FR[pos0+1] = (normal_x*Dm[1]-normal_y*Dm[2]);
        FR[pos0+2] = (normal_y*Dm[1]+normal_x*Dm[2]);
        FR[pos0+3] = Dm[3];
        
        FL[pos1]   = Dp[0];
        FL[pos1+1] = (normal_x*Dp[1]-normal_y*Dp[2]);
        FL[pos1+2] = (normal_y*Dp[1]+normal_x*Dp[2]);
        FL[pos1+3] = Dp[3];
    }
    else{
        if(pos1==-1){
            FL[pos0]   = Dm[0];
            FL[pos0+1] = (normal_x*Dm[1]-normal_y*Dm[2]);
            FL[pos0+2] = (normal_y*Dm[1]+normal_x*Dm[2]);
            FL[pos0+3] = Dm[3];           
        }
        else{
            FR[pos0]   = Dm[0];
            FR[pos0+1] = (normal_x*Dm[1]-normal_y*Dm[2]);
            FR[pos0+2] = (normal_y*Dm[1]+normal_x*Dm[2]);
            FR[pos0+3] = Dm[3];
        }

    }
    
#endif
 
    
    
}

#endif

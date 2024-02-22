#ifndef _MATRIZ_H_
#define _MATRIZ_H_

#include "Constantes.hxx"


__device__ double pow4(double x) {
    return x*x*x*x;
}

__device__ void calculo_u(double h, double2 q, double epsilon_h, double2& u) {

    double modu;
    u.x = q.x/(h+EPSILON)*(h>=EPS);
    u.y = q.y/(h+EPSILON)*(h>=EPS);
    modu = sqrt(u.x*u.x+u.y*u.y);
    if (modu>VMAX) {
        u.x=u.x/modu*VMAX;
        u.y=u.y/modu*VMAX;
    }
}






__device__ void filtro_estado(double& h, double2& q , double epsilon_h)
{
    double2 u;
    
    h = h*(h>0.0);
    
    calculo_u(h,q,epsilon_h,u);
    
    q.x = h*u.x;
    q.y = h*u.y;
}


__device__ void evalua_rh(double h,double dh,double val, double& hr) {    
    hr=h+val*dh;
}

__device__ void evalua_rzb(double zc,double dz,double val, double& z) {    
    z=zc+val*dz;
}

__device__ void evalua_rq(double2 q,double2 dq, double val, double2& qr) {
    qr.x=q.x+val*dq.x;
    qr.y=q.y+val*dq.y;
    
}


__device__ double bedload(double h, double2 un, double db, double dx){
	double qb;
	double K1;
	double O;
	double u;
	double aux;
    
    u = un.x;
    double modu = sqrt(un.x*un.x+un.y*un.y) ;


// #if bedload_formula==0

    O = d_Opar*modu*modu;

    K1 = 8.0*u/(modu+EPSILON);
    aux=max(O-d_OC,0.0);
    qb=K1*d_Q*pow(aux,1.5)*(h>5*d_DS)*(h>d_HQB);

    // Ashida & Michiue
    K1 = 17.0*u/(modu+EPSILON);
    qb=K1*d_Q*max(O-d_OC,0.0)*max(sqrt(O)-sqrt(d_OC),0.0)*(h>d_HQB)*(h>5*d_DS);
    
    // qb *= 53;
// #else
//
//     double epsb, epss, W, phi;
//     double w, tanphi, tanbeta;
//
//     epsb=0.13;
//     epss=0.01;
//     W = 0.04;
//     w = 1000*5e-3*abs(modu*modu*modu);
//     tanphi = tan(25*PI/180);
//     tanbeta= db/dx;
//
//     qb = 0.0001*(epsb/(tanphi-tanbeta) + epss/((W/modu)-tanbeta))*w;
//
//     // if(abs(qb)>1e-5)
//     // printf("%.15f\n",qb);
//
// #endif

    
    return -qb;  // La explicacion de este signo se encuentra en apunte 1 del resolvedor de Riemann.
}


__device__ void flujoC(double h, double2 u, double *F,int fixed_sediment, double db, double dx){
    
    F[0] = h*u.x;
    F[1] = F[0]*u.x;
    F[2] = F[0]*u.y;
    F[3] = bedload(h,u, db, dx) *fixed_sediment  * SEDIMENT;

}

__device__ void roe_average(double h0, double2 u0, double z0, double h1, double2 u1, double z1, double &hm, double2 &um, double &zm){
    
    hm = 0.5*(h0+h1);
    
    um.x = (sqrt(h0)*u0.x + sqrt(h1)*u1.x)/(sqrt(h0)+sqrt(h1));
    um.y = (sqrt(h0)*u0.y + sqrt(h1)*u1.y)/(sqrt(h0)+sqrt(h1));
    
    zm = 0.5*(z0+z1);
}


__device__ void calc_autoval_aprox0(double h, double2 u, double *autovalores){
	double c=sqrt(G*h);
    
    autovalores[0] = u.x-c;
    autovalores[1] = 0.0;
    autovalores[2] = u.x+c;

}


#endif

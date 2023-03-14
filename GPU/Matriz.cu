#ifndef _MATRIZ_H_
#define _MATRIZ_H_

#include "Constantes.hxx"

// extern __device__ double d_OC;
// extern __device__ double d_Q;
// extern __device__ double d_DS;
// extern __device__ double d_Opar;
// extern __device__ double d_HQB;


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


__device__ double bedload(double h, double2 un){
	double qb;
	double K1;
	double O;
	double u;
	double aux;
    
    u = un.x;
    double modu = sqrt(un.x*un.x+un.y*un.y) ;
	
// 	u=div(U[1],U[0]);
// 	O=d_Opar*u*u;
    O = d_Opar*modu*modu;

#if model == 0
//     K1=8.0*SGN(u);
    K1 = 8.0*u/(modu+EPSILON);
    aux=max(O-d_OC,0.0);
    qb=K1*d_Q*pow(aux,1.5)*(h>5*d_DS)*(h>d_HQB);
    
#elif model==1
//     K1=17.0*SGN(u);
    K1 = 17.0*u/modu;
    qb=K1*d_Q*max(O-d_OC,0.0)*max(sqrt(O)-sqrt(d_OC),0.0)*(h>d_HQB);
#endif

    
    return -qb ;  // La explicacion de este signo se encuentra en apunte 1 del resolvedor de Riemann.
}


__device__ void flujoC(double h, double2 u, double *F){
    
    F[0] = h*u.x;
    F[1] = F[0]*u.x;
    F[2] = F[0]*u.y;
    F[3] = bedload(h,u);
    
}

__device__ void roe_average(double h0, double2 u0, double z0, double h1, double2 u1, double z1, double &hm, double2 &um, double &zm){
    
    hm = 0.5*(h0+h1);
    
    um.x = (sqrt(h0)*u0.x + sqrt(h1)*u1.x)/(sqrt(h0)+sqrt(h1));
    um.y = (sqrt(h0)*u0.y + sqrt(h1)*u1.y)/(sqrt(h0)+sqrt(h1));
    
//     qm.x = 9.0;
//     qm.y = 3.0;
//     qm.x *= hm;
//     qm.y *= hm;

    zm = 0.5*(z0+z1);
}


__device__ void calc_autoval_aprox0(double h, double2 u, double *autovalores){
	double c=sqrt(G*h);

//     Entiendo que aqui tengo que tomar el modulo de u. Coger el mas restrictivo podria ser excesivo en problemas donde una direccion predomine sobre otra.
    
//     double um = sqrt(u.x*u.x+u.y*u.y);
//     um = u.x/um;
//     autovalores[0]=um-c;
// 	autovalores[1]=0.0;
// 	autovalores[2]=um+c;
    
    autovalores[0] = u.x-c;
    autovalores[1] = 0.0;
    autovalores[2] = u.x+c;
    
//     autovalores[0] = min(u.x-c, u.y-c);
//     autovalores[1] = 0;
//     autovalores[2] = max(u.x+c, u.y+c);
// 	
}


#endif

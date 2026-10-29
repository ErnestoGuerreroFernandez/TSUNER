// TSUNER Copyright (C) 2025 Ernesto Guerrero Fernandez
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.



#pragma once
#include "GPU_kernels.cuh"


// static __device__ __forceinline__ inline double pow4(double x) {
//     return x*x*x*x;
// }

__device__ __forceinline__ void u_calculator(double h, double2 q, double epsilon_h, double2& u) {

    double modu;
    u.x = q.x/(h+EPSILON)*(h>=EPS);
    u.y = q.y/(h+EPSILON)*(h>=EPS);
    modu = sqrt(u.x*u.x+u.y*u.y);
    if (modu>VMAX) {
        u.x=u.x/modu*VMAX;
        u.y=u.y/modu*VMAX;
    }
}


__device__ __forceinline__ void u_calculator_double(double h, double qx, double qy, double epsilon_h, double2& u) {

    double modu;
    u.x = qx/(h+EPSILON)*(h>=EPS);
    u.y = qy/(h+EPSILON)*(h>=EPS);
    modu = sqrt(u.x*u.x+u.y*u.y);
    if (modu>VMAX) {
        u.x=u.x/modu*VMAX;
        u.y=u.y/modu*VMAX;
    }
}

__device__ __forceinline__ void calculate_erosion(double2 u, double modu, double Rp, double vs, const SimParams &params, double &Fe){

    double Es, Z, Z5;
    double alpha1, alpha2;

    if (Rp>2.36){
        alpha1 = 1;
        alpha2 = 0.6;
    }
    else{
        alpha1 = 0.586;
        alpha2 = 1.23;
    }

    Z = alpha1 * pow(Rp,alpha2) * (modu*sqrt(params.FR))/vs;
    Z5 = Z*Z*Z*Z*Z;

    Es = (1.3e-7*Z5)/(1+4.3e-7*Z5);

    Fe = vs * (1-params.porosity) *Es;

    if(Fe>0.1){
        Fe = 0.1;
    }
    Fe = Fe*(Fe>0.0);

}

__device__ __forceinline__ void calculate_deposition(double c, double vs, double Rp, const SimParams &params, double &Fd){

    double Cb;

    Cb = (1.64+0.4*pow(params.DS/params.Dsg,1.64))*c;

    // double m;
    // m = 4.45 * pow(Rp,-0.1);
    // Fd = vs * Cb * pow((1-Cb),m);
    Fd = vs*Cb;
    if(Fd>0.1){
        Fd = 0.1;
    }
    Fd = Fd*(Fd>0.0);

    // printf("%.14f  %.14f \n ",Fd,c);




}



__device__ __forceinline__ void state_filter(double& h, double2& q, double &hc, double epsilon_h)
{
    double2 u;
    double c;
    
    h = h*(h>0.0);
    c = hc/(h+EPSILON);
    c = c*(c>0.0);

    u_calculator(h,q,epsilon_h,u);
    
    q.x = h*u.x;
    q.y = h*u.y;
    hc = h*c;
}
__device__ __forceinline__ void state_filter_double(double& h, double& qx, double& qy, double &hc, double epsilon_h)
{
    double2 u;
    double c;

    h = h*(h>0.0);
    c = hc/(h+EPSILON);
    c = c*(c>0.0);

    u_calculator_double(h,qx,qy,epsilon_h,u);

    qx = h*u.x;
    qy = h*u.y;
    hc = h*c;
}


__device__ __forceinline__ void rh_eval(double h,double dh,double val, double& hr) {
    hr=h+val*dh;
}

__device__ __forceinline__ void rzb_eval(double zc,double dz,double val, double& z) {
    z=zc+val*dz;
}

__device__ __forceinline__ void rq_eval(double2 q,double2 dq, double val, double2& qr) {
    qr.x=q.x+val*dq.x;
    qr.y=q.y+val*dq.y;
    
}


__device__ __forceinline__ double bedload(double h, double2 un, double db, const SimParams &params){
	double qb;
	double K1;
	double O;
	double u;
    
    u = un.x;
    double modu = sqrt(un.x*un.x+un.y*un.y) ;



    O = params.Opar*modu*modu;

    // Ashida & Michiue
    K1 = 17.0*u/(modu+EPSILON);
    qb=K1*params.Q*max(O-params.OC,0.0)*max(sqrt(O)-sqrt(params.OC),0.0)*(h>params.HQB)*(h>5*params.DS);
    
    return -qb;
}


__device__ __forceinline__ void Cflux(double h, double2 u, double hc, double *F, int fixed_sediment, double db, const SimParams &params){
    
    F[0] = h*u.x;
    F[1] = F[0]*u.x;
    F[2] = F[0]*u.y;
    F[3] = bedload(h,u, db, params) *fixed_sediment  * BEDLOAD;
    F[4] = hc*u.x;

}

__device__ __forceinline__ void roe_average(double h0, double2 u0, double z0, double h1, double2 u1, double z1, double &hm, double2 &um, double &zm){
    
    hm = 0.5*(h0+h1);
    
    um.x = (sqrt(h0)*u0.x + sqrt(h1)*u1.x)/(sqrt(h0)+sqrt(h1));
    um.y = (sqrt(h0)*u0.y + sqrt(h1)*u1.y)/(sqrt(h0)+sqrt(h1));
    
    zm = 0.5*(z0+z1);
}


__device__ __forceinline__ void calc_eigen_aprox0(double h, double2 u, double *autovalores){
	double c=sqrt(G*h);
    
    autovalores[0] = u.x-c;
    autovalores[1] = 0.0;
    autovalores[2] = u.x+c;

}


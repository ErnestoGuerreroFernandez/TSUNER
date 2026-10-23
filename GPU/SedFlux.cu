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



#ifndef _MATRIZ_H_
#define _MATRIZ_H_

#include "Sim_parameters.hxx"


__device__ double pow4(double x) {
    return x*x*x*x;
}

__device__ void u_calculator(double h, double2 q, double epsilon_h, double2& u) {

    double modu;
    u.x = q.x/(h+EPSILON)*(h>=EPS);
    u.y = q.y/(h+EPSILON)*(h>=EPS);
    modu = sqrt(u.x*u.x+u.y*u.y);
    if (modu>VMAX) {
        u.x=u.x/modu*VMAX;
        u.y=u.y/modu*VMAX;
    }
}






__device__ void state_filter(double& h, double2& q , double epsilon_h)
{
    double2 u;
    
    h = h*(h>0.0);
    
    u_calculator(h,q,epsilon_h,u);
    
    q.x = h*u.x;
    q.y = h*u.y;
}


__device__ void rh_eval(double h,double dh,double val, double& hr) {
    hr=h+val*dh;
}

__device__ void rzb_eval(double zc,double dz,double val, double& z) {    
    z=zc+val*dz;
}

__device__ void rq_eval(double2 q,double2 dq, double val, double2& qr) {
    qr.x=q.x+val*dq.x;
    qr.y=q.y+val*dq.y;
    
}


__device__ double bedload(double h, double2 un, double db){
	double qb;
	double K1;
	double O;
	double u;
    
    u = un.x;
    double modu = sqrt(un.x*un.x+un.y*un.y) ;



    O = d_Opar*modu*modu;

    // Ashida & Michiue
    K1 = 17.0*u/(modu+EPSILON);
    qb=K1*d_Q*max(O-d_OC,0.0)*max(sqrt(O)-sqrt(d_OC),0.0)*(h>d_HQB)*(h>5*d_DS);
    
    return -qb;
}


__device__ void Cflux(double h, double2 u, double *F,int fixed_sediment, double db){
    
    F[0] = h*u.x;
    F[1] = F[0]*u.x;
    F[2] = F[0]*u.y;
    F[3] = bedload(h,u, db) *fixed_sediment  * SEDIMENT;

}

__device__ void roe_average(double h0, double2 u0, double z0, double h1, double2 u1, double z1, double &hm, double2 &um, double &zm){
    
    hm = 0.5*(h0+h1);
    
    um.x = (sqrt(h0)*u0.x + sqrt(h1)*u1.x)/(sqrt(h0)+sqrt(h1));
    um.y = (sqrt(h0)*u0.y + sqrt(h1)*u1.y)/(sqrt(h0)+sqrt(h1));
    
    zm = 0.5*(z0+z1);
}


__device__ void calc_eigen_aprox0(double h, double2 u, double *autovalores){
	double c=sqrt(G*h);
    
    autovalores[0] = u.x-c;
    autovalores[1] = 0.0;
    autovalores[2] = u.x+c;

}


#endif

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



#include "GPU_kernels.cuh"
#include "SedFlux.cuh"

static __device__ inline double avg(double a, double b)
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

// static __device__ inline double minmod(double a, double b, double c)
// {
//     double d;
//     if ((a>0.0) && (b>0.0) && (c>0.0)) {
//         d=min(a,min(b,c));
//     }
//     else {
//         if ((a<0.0) && (b<0.0) && (c<0.0)) {
//             d=max(a,max(b,c));
//         }
//         else {
//             d=0.0;
//         }
//     }
//     return d;
// }

static __device__ inline void friction(double* d_Vk, double* d_V0, int pos, double epsilon_h, double bottom_friction, double dt) {
    
    double h;
    double2 q, u;
    
    h = d_V0[pos];
    q.x = d_V0[pos+1];
    q.y = d_V0[pos+2];
    
    u_calculator(h,q,epsilon_h,u);
    double modu = sqrt(u.x*u.x+u.y*u.y);
    
    double c = dt*bottom_friction*modu/(h+EPSILON);
    
#if Manning ==1
    c = dt*bottom_friction*bottom_friction*modu/(pow(h,4.0/3.0));
    d_Vk[pos+1] *= c;
    d_Vk[pos+2] *= c;
#else
    d_Vk[pos+1] /= (1.0+c);
    d_Vk[pos+2] /= (1.0+c);
#endif    
    
            
}

__device__ void erosion_deposition(double* d_Vk, double* d_V0, const SimParams &params, double dt, int pos){

    double Fe, Fd, phi;
    double h, c;
    double vs, Rp, nu;
    double2 q, u;

    nu = 1.002e-6; // Water kinematic viscosity
    double s = params.R-1; // Diff of relative densities

    h = d_V0[pos];
    q.x = d_V0[pos+1];
    q.y = d_V0[pos+2];
    c = d_V0[pos+4]/(h+params.epsilon_h);



    u_calculator(h,q,params.epsilon_h,u);
    double modu = sqrt(u.x*u.x+u.y*u.y);

    vs = sqrt( (13.95*nu/params.DS)*(13.95*nu/params.DS) + 1.09 * s * G* params.DS ) - 13.95*nu/params.DS  ;
    Rp = (params.DS*sqrt(s*G*params.DS))/nu; // Particle Reynolds number

    calculate_erosion(u,modu,Rp,vs,params,Fe);
    calculate_deposition(c,vs,params,Fd);

    phi = (Fe - Fd) * SUSPENSION;


    // double aux = Vhc;
    // double aux2 = Vb;
    // Vhc = Vhc + dt*phi;//*(h>0.001);
    // Vb = Vb - dt*phi;//*(h>0.001);

    d_Vk[pos+3] = d_Vk[pos+3] - dt*phi;
    d_Vk[pos+4] = d_Vk[pos+4] + dt*phi;



}




__global__ void reconstructor(double* d_V, double* reconsL, double* reconsR, double* reconsD, double* reconsU, double delta_T, double inv_dx, double inv_dy, SimParams params, double* d_Vk, double *d_V0) {


    int i,j;
    int pos;
	i = blockIdx.x*NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
	j = blockIdx.y*NUM_THREADS_VERTICAL_VOL + threadIdx.y;

    if ((i<params.npx) && (j<params.npy)) {
        pos = (j*params.npx + i)*5;
        
#if order == 1


        reconsL[pos]   = reconsR[pos]   = reconsU[pos]   = reconsD[pos]   = d_V[pos];
        reconsL[pos+1] = reconsR[pos+1] = reconsU[pos+1] = reconsD[pos+1] = d_V[pos+1];
        reconsL[pos+2] = reconsR[pos+2] = reconsU[pos+2] = reconsD[pos+2] = d_V[pos+2];
        reconsL[pos+3] = reconsR[pos+3] = reconsU[pos+3] = reconsD[pos+3] = d_V[pos+3];
        reconsL[pos+4] = reconsR[pos+4] = reconsU[pos+4] = reconsD[pos+4] = d_V[pos+4];



        
#else   

        int posleft,posright;
        int im,ip,jm,jp;
        double dta,Hl,Hi,Hr;
        double hi,hl,hr;
        double2 ui,ul,ur;
        double dl, dr, du;
        double ci, cr, cl;
        double deta;
        
        // **********************
        // x direction
        // **********************

        hi=d_V[pos];
        Hi=d_V[pos+3];
        ci=d_V[pos+4]/(hi+EPSILON);

        im=i-1*(i>0);
        ip=i+1*(i<params.npx-1);
#if WBC==2
        if(i==0) im = params.npx-1;
        if(i==params.npx-1) ip = 0;
#endif
        posleft=(j*params.npx + im)*5;
        posright=(j*params.npx + ip)*5;
        hl=d_V[posleft];
        Hl=d_V[posleft+3];
        cl=d_V[posleft+4]/(hl+EPSILON);
        hr=d_V[posright];
        Hr=d_V[posright+3];
        cr=d_V[posright+4]/(hr+EPSILON);
        dta=delta_T*inv_dx;

        
        // h reconstruction
        dl = hi - hl;
        dr = hr - hi;
        double dh = avg(dl,dr);
        if(0.5*fabs(dh)>hi){
            dh = 2.0*hi*SGN(dh);
        }
        reconsL[pos] = hi - 0.5*dh;
        reconsR[pos] = hi + 0.5*dh;
        
        // eta and zb reconstruction
        double zm;
        zm = min(Hl,Hi);
        dl = max(hi-Hi+zm,0.0) - max(hl-Hl+zm,0.0);
        zm = min(Hr,Hi);
        dr = max(hr-Hr+zm,0.0) - max(hi-Hi+zm,0.0);
        deta = avg(dl,dr);
        double dz = dh - deta;
        
        reconsL[pos+3] = Hi - 0.5*dz;
        reconsR[pos+3] = Hi + 0.5*dz;
        
        // u reconstruction
        u_calculator_double(hi,d_V[pos+1],d_V[pos+2], params.epsilon_h,ui);
        u_calculator_double(hl,d_V[posleft+1],d_V[posleft+2],params.epsilon_h,ul);
        u_calculator_double(hr,d_V[posright+1],d_V[posright+2],params.epsilon_h,ur);
        dl = ui.x - ul.x;
        dr = ur.x - ui.x;
        du = avg(dl,dr);
        double dq = hi*du + ui.x * dh;


        
        reconsL[pos+1] = d_V[pos+1] - 0.5*dq;
        reconsR[pos+1] = d_V[pos+1] + 0.5*dq;
        
        
        dl = ui.y - ul.y;
        dr = ur.y - ui.y;
        du = avg(dl,dr);
        dq = hi*du + ui.y * dh;
        
        reconsL[pos+2] = d_V[pos+2] - 0.5*dq;
        reconsR[pos+2] = d_V[pos+2] + 0.5*dq;


        // hc reconstruction
        dl = ci - cl;
        dr = cr - ci;
        du = avg(dl,dr);

        reconsL[pos+4] = (ci - 0.5*du)*reconsL[pos];
        reconsR[pos+4] = (ci + 0.5*du)*reconsR[pos];
        
        
        // Integral
        d_Vk[pos+1] -= dta*G*hi*deta;
        
#if RK==0
        double Tpx = G*hi*deta;
#endif
        
        // **********************
        // y direction
        // **********************
        

        jm=j-1*(j>0);
        jp=j+1*(j<params.npy-1);
#if WBC==2
        if(j==0) jm = params.npy-1;
        if(j==params.npy-1) jp = 0;
#endif
        posleft=(jm*params.npx +i)*5;
        posright=(jp*params.npx +i)*5;
        hl=d_V[posleft];
        Hl=d_V[posleft+3];
        hr=d_V[posright];
        Hr=d_V[posright+3];
        dta=delta_T*inv_dy;
        
        // h reconstruction
        dl = hi - hl;
        dr = hr - hi;
        dh = avg(dl,dr);
        if(0.5*fabs(dh)>hi){
            dh = 2.0*hi*SGN(dh);
        }
        reconsD[pos] = hi - 0.5*dh;
        reconsU[pos] = hi + 0.5*dh;
        

        
        // eta and zb reconstruction
        zm = min(Hl,Hi);
        dl = max(hi-Hi+zm,0.0) - max(hl-Hl+zm,0.0);
        zm = min(Hr,Hi);
        dr = max(hr-Hr+zm,0.0) - max(hi-Hi+zm,0.0);
        deta = avg(dl,dr);
        dz = dh - deta;
        
        reconsD[pos+3] = Hi - 0.5*dz;
        reconsU[pos+3] = Hi + 0.5*dz;
        

        // u reconstruction
        u_calculator_double(hl,d_V[posleft+1],d_V[posleft+2],params.epsilon_h,ul);
        u_calculator_double(hr,d_V[posright+1],d_V[posright+2],params.epsilon_h,ur);
        dl = ui.x - ul.x;
        dr = ur.x - ui.x;
        du = avg(dl,dr);
        dq = hi*du + ui.x * dh;
        
        reconsD[pos+1]  = d_V[pos+1] - 0.5*dq;
        reconsU[pos+1]  = d_V[pos+1] + 0.5*dq;
        
        dl = ui.y - ul.y;
        dr = ur.y - ui.y;
        du = avg(dl,dr);
        dq = hi*du + ui.y * dh;
        
        reconsD[pos+2]  = d_V[pos+2] - 0.5*dq;
        reconsU[pos+2]  = d_V[pos+2] + 0.5*dq;

        // hc reconstruction
        dl = ci-cl;
        dr = cr-ci;
        du = avg(dl,dr);

        reconsD[pos+4] = (ci-0.5*du)*reconsD[pos];
        reconsU[pos+4] = (ci+0.5*du)*reconsU[pos];

        
        // Integral
        d_Vk[pos+2] -= dta*G*hi*deta;
        
#if RK==0
        double Tpy = G*hi*deta;
#endif
        
                
#if RK==0            
        double2 uRR, uL, uU, uD;
        u_calculator_double(reconsR[pos],reconsR[pos+1],reconsR[pos+2],params.epsilon_h,uRR);
        u_calculator_double(reconsL[pos],reconsL[pos+1],reconsL[pos+2],params.epsilon_h,uL);
        u_calculator_double(reconsU[pos],reconsU[pos+1],reconsU[pos+2],params.epsilon_h,uU);
        u_calculator_double(reconsD[pos],reconsD[pos+1],reconsD[pos+2],params.epsilon_h,uD);
    
        double qhx, qhy;
        qhx = (reconsR[pos+1]-reconsL[pos+1]) *inv_dx;
        qhy = (reconsU[pos+2]-reconsD[pos+2]) *inv_dy;
        
        
        // MUSCL-Hancock
        reconsL[pos] -= 0.5 * delta_T * (qhx + qhy);
        reconsR[pos] -= 0.5 * delta_T * (qhx + qhy);
        reconsD[pos] -= 0.5 * delta_T * (qhx + qhy);
        reconsU[pos] -= 0.5 * delta_T * (qhx + qhy);
        

        qhx = (reconsR[pos+1]*uRR.x-reconsL[pos+1]*uL.x)*inv_dx;
        qhy = (reconsU[pos+1]*uU.y-reconsD[pos+1]*uD.y)*inv_dy;
        reconsL[pos+1] -= 0.5 * delta_T * (qhx + qhy);
        reconsR[pos+1] -= 0.5 * delta_T * (qhx + qhy);
        reconsD[pos+1] -= 0.5 * delta_T * (qhx + qhy);
        reconsU[pos+1] -= 0.5 * delta_T * (qhx + qhy);
        

        qhx = (reconsR[pos+1]*uRR.y - reconsL[pos+1]*uL.y) *inv_dx;
        qhy = (reconsU[pos+2]*uU.y - reconsD[pos+2]*uD.y) *inv_dy;
        reconsL[pos+2] -= 0.5 * delta_T * (qhx + qhy);
        reconsR[pos+2] -= 0.5 * delta_T * (qhx + qhy);
        reconsD[pos+2] -= 0.5 * delta_T * (qhx + qhy);
        reconsU[pos+2] -= 0.5 * delta_T * (qhx + qhy);
        

        qhx = (bedload(reconsR[pos],uRR,0,params)-bedload(reconsL[pos],uL,0,params))*inv_dx;
        qhy = (bedload(reconsU[pos],uU,0,params)-bedload(reconsD[pos],uD,0,params))*inv_dy;
        reconsL[pos+3] -= 0.5 * delta_T * (qhx + qhy);
        reconsR[pos+3] -= 0.5 * delta_T * (qhx + qhy);
        reconsD[pos+3] -= 0.5 * delta_T * (qhx + qhy);
        reconsU[pos+3] -= 0.5 * delta_T * (qhx + qhy);


        qhx = reconsR[pos+4]*uRR.x - reconsL[pos+4]*uL.x;
        qhy = reconsU[pos+4]*uU.y - reconsD[pos+4]*uD.y;
        reconsL[pos+4] -= 0.5 * delta_T * (qhx + qhy);
        reconsR[pos+4] -= 0.5 * delta_T * (qhx + qhy);
        reconsD[pos+4] -= 0.5 * delta_T * (qhx + qhy);
        reconsU[pos+4] -= 0.5 * delta_T * (qhx + qhy);
        
        
        
        reconsL[pos+1] -= 0.5 * delta_T*inv_dx * Tpx;
        reconsU[pos+1] -= 0.5 * delta_T*inv_dx * Tpx;
        reconsD[pos+1] -= 0.5 * delta_T*inv_dx * Tpx;
        reconsR[pos+1] -= 0.5 * delta_T*inv_dx * Tpx;
        
        reconsL[pos+2] -= 0.5 * delta_T*inv_dy * Tpy;
        reconsU[pos+2] -= 0.5 * delta_T*inv_dy * Tpy;
        reconsD[pos+2] -= 0.5 * delta_T*inv_dy * Tpy;
        reconsR[pos+2] -= 0.5 * delta_T*inv_dy * Tpy;

#endif   
        
        
        // Filters
        state_filter_double(reconsL[pos],reconsL[pos+1],reconsL[pos+2],reconsL[pos+4],params.epsilon_h);
        state_filter_double(reconsR[pos],reconsR[pos+1],reconsR[pos+2],reconsR[pos+4],params.epsilon_h);
        state_filter_double(reconsU[pos],reconsU[pos+1],reconsU[pos+2],reconsU[pos+4],params.epsilon_h);
        state_filter_double(reconsD[pos],reconsD[pos+1],reconsD[pos+2],reconsD[pos+4],params.epsilon_h);
        
        
#endif
        
        
        // Frictions
        friction(d_Vk, d_V0, pos, params.epsilon_h, params.FR, delta_T);
        erosion_deposition(d_Vk, d_V0, params, delta_T, pos);

    }

}

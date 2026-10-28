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


__global__ void GetDeltaT(double* d_V, double* d_deltaTV, SimParams params)
{
    int pos,posv,i,j;
    double2 q;
    double2 u;
    double h;
    double cg,lx,ly,dt;
    

    i = blockIdx.x*NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
	j = blockIdx.y*NUM_THREADS_VERTICAL_VOL + threadIdx.y;
    

    if ((i<params.npx) && (j<params.npy)) {
        
        pos=j*params.npx+i;
        posv = pos*5;
        h=d_V[posv];
        q.x=d_V[posv+1];
        q.y=d_V[posv+2];
        u_calculator(h, q, params.epsilon_h, u);
        
        cg=sqrt(G*h);
        lx=fabs(u.x)+cg;
        ly=fabs(u.y)+cg;
//         dt=CFL*0.5*dx*dy/(dy*lx+dx*ly+EPSILON);
        dt = params.CFL * 0.5 * min(params.dx/lx, params.dy/ly);
        
        d_deltaTV[pos]=dt;



    }
}

static __device__ inline void interpolate_timeSPL(double *CC, double *tCC, long unsigned int timeCC, double current_time, double &h1, double &zb1, double2 &q1, int npx, int npy, int pos){
    
    
    double dt;
    dt = tCC[1]-tCC[0];
    int t, tp;
    t = int((current_time-tCC[0])/dt);
    t*=(t>0);
    if(t>timeCC-1)
        t=timeCC-1;
    tp = t+1;
    int posCC = t*(2*npx+2*npy)*4 + pos*4;
    int posCCp = tp*(2*npx+2*npy)*4 + pos*4;
    
    double m, b;
    m = (CC[posCCp]-CC[posCC])/(tCC[tp]-tCC[t]);
    b = CC[posCC] - m*tCC[t];
    h1 = m*current_time+b;
    h1*=UNITS;
    
    m = (CC[posCCp+3]-CC[posCC+3])/(tCC[tp]-tCC[t]);
    b = CC[posCC+3] - m*tCC[t];
    zb1 = m*current_time+b;

    h1 += zb1;    

    m = (CC[posCCp+1]-CC[posCC+1])/(tCC[tp]-tCC[t]);
    b = CC[posCC+1] - m*tCC[t];
    q1.x = m*current_time+b;    
    q1.x*= UNITS* h1;

    m = (CC[posCCp+2]-CC[posCC+2])/(tCC[tp]-tCC[t]);
    b = CC[posCC+2] - m*tCC[t];
    q1.y = m*current_time+b;   
    q1.y*= UNITS*h1;

    
}





__global__ void hmax_calculator(double* d_hmax, double * d_V, double epsilon_h, int npx, int npy) {
    int i,j;
    int pos,posv;
    i = blockIdx.x*NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    j = blockIdx.y*NUM_THREADS_VERTICAL_VOL + threadIdx.y;
    if ((i<npx) && (j<npy)) {
       pos=j*npx+i;
       posv = pos*5;
       if (d_V[posv]>epsilon_h) {
         d_hmax[pos]=max(d_hmax[pos],d_V[posv]);
       } 
    }

}
__global__ void bmin_calculator(double* d_bmin, double * d_V, int npx, int npy) {
    int i,j;
    int pos,posv;
    i = blockIdx.x*NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    j = blockIdx.y*NUM_THREADS_VERTICAL_VOL + threadIdx.y;
    if ((i<npx) && (j<npy)) {
        pos=j*npx+i;
        posv=pos*5;
        d_bmin[pos] = max(d_bmin[pos],d_V[posv+3]);
    }
}

__global__ void eta_max_calculator(double* d_eta_max, double *d_V, double epsilon_h, int npx, int npy) {
    int i, j;
    int pos,posv;
    i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    j = blockIdx.y * NUM_THREADS_VERTICAL_VOL + threadIdx.y;
    if ((i < npx) && (j < npy)) {
        pos = j * npx + i;
        posv = pos*5;
        if (d_V[posv] > epsilon_h) {
            d_eta_max[pos] = max(d_eta_max[pos], d_V[posv] - d_V[posv+3]);
        }
    }
}

__global__ void set_inicializacion(double* d_hmax, int npx, int npy) {
    int i,j;
    int pos;
    i = blockIdx.x*NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    j = blockIdx.y*NUM_THREADS_VERTICAL_VOL + threadIdx.y;
    if ((i<npx) && (j<npy)) {
       pos=j*npx+i;
       d_hmax[pos] = -1e30;
    }
}
__global__ void set_inicializacion_bmin(double *d_bmin, int npx, int npy){
    int i,j;
    int pos;
    i = blockIdx.x*NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    j = blockIdx.y*NUM_THREADS_VERTICAL_VOL + threadIdx.y;
    if ((i<npx) && (j<npy)) {
        pos=j*npx+i;
        d_bmin[pos] = -1e30;
    }
}


__global__ void set_inicializacion_eta_max(double *d_eta_max, int npx, int npy) {
    int i, j;
    int pos;
    i = blockIdx.x * NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
    j = blockIdx.y * NUM_THREADS_VERTICAL_VOL + threadIdx.y;
    if ((i < npx) && (j < npy)) {
        pos = j * npx + i;
        d_eta_max[pos] = -1e30;
    }
}



__device__ void spongeLayer(double *d_Ah, double2 *d_Ahu, double *d_Azb, int npx, int npy, double *d, int pos, int pos1, int pos2, int posd1, int posd2, double *CC, double *tCC, long unsigned int timeCC, double current_time){
    
    double h1, h2, H1, H2;
    double2 q1, q2;
    
        
    // Corner
    if(pos2>0){
        interpolate_timeSPL(CC, tCC, timeCC, current_time, h1, H1, q1, npx, npy, pos1);
        interpolate_timeSPL(CC, tCC, timeCC, current_time, h2, H2, q2, npx, npy, pos2);
        double aux = 0.5*(d[posd1]+d[posd2])*0.5;
        d_Ah[pos] = aux*(h1+h2) + (1.0-d[posd1])*d_Ah[pos]*0.5 + (1.0-d[posd2])*d_Ah[pos] *0.5 ;
        d_Ahu[pos].x = aux*(q1.x+q2.x) + (1.0-d[posd1])*d_Ahu[pos].x*0.5 + (1.0-d[posd2])*d_Ahu[pos].x *0.5 ;
        d_Ahu[pos].y = aux*(q1.y+q2.y) + (1.0-d[posd1])*d_Ahu[pos].y*0.5 + (1.0-d[posd2])*d_Ahu[pos].y *0.5 ;
        d_Azb[pos] = aux*(H1+H2) + (1.0-d[posd1])*d_Azb[pos]*0.5 + (1.0-d[posd2])*d_Azb[pos] *0.5 ;
    }

    else{
        interpolate_timeSPL(CC, tCC, timeCC, current_time, h1, H1, q1, npx, npy, pos1);
        d_Ah[pos]    = d[posd1]*h1 + (1.0-d[posd1])*d_Ah[pos];
        d_Ahu[pos].x = d[posd1]*q1.x + (1.0-d[posd1])*d_Ahu[pos].x;
        d_Ahu[pos].y = d[posd1]*q1.y + (1.0-d[posd1])*d_Ahu[pos].y;
        d_Azb[pos]   = d[posd1]*H1 + (1.0-d[posd1])*d_Azb[pos];
    }


}


// Untested!
__device__ double dynamicRuptureInterpolation(double *dyrup, float *time_dyrup, double current_time, int pos, int npx, int npy){

    // return 0.0;
    float dt;
    dt = time_dyrup[1]-time_dyrup[0];
    int t, tp;
    t = int((current_time-time_dyrup[0])/dt);
    // printf("%d  %f \n",t,current_time);
    t*=(t>0);
    if(t>40-1)
        return 0.0;
    tp = t+1;
    int posCC = t*npx*npy + pos;
    int posCCp = tp*npx*npy + pos;

    double m, b;
    m = (dyrup[posCCp]-dyrup[posCC])/(time_dyrup[tp]-time_dyrup[t]);
    b = dyrup[posCC] - m*time_dyrup[t];
    return (m*current_time+b);

}




__global__ void states_updater(double* d_Vk, double* d_V0, double* d_V, int npx, int npy, int k, double epsilon_h, double invx, double invy, double *FL, double *FR, double *FU, double *FD, double *d, double current_time, double *d_CC, double *d_tCC, long unsigned int timeCC, double *dyrup, float *time_dyrup)
{
    int i,j;
    int pos;
    double coef0,coef1;
    double2 q;
    double h, z, hc;

    i = blockIdx.x*NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
	j = blockIdx.y*NUM_THREADS_VERTICAL_VOL + threadIdx.y;

    if ((i<npx) && (j<npy)) {
        
        pos = j*npx + i;
        int posF = pos*5;


        
        d_Vk[posF]   -= invx*(FL[posF]  +FR[posF])    + invy*(FD[posF]+  FU[posF]);
        d_Vk[posF+1] -= invx*(FL[posF+1]+FR[posF+1])  + invy*(FD[posF+1]+FU[posF+1]);
        d_Vk[posF+2] -= invx*(FL[posF+2]+FR[posF+2])  + invy*(FD[posF+2]+FU[posF+2]);
        d_Vk[posF+3] -= invx*(FL[posF+3]+FR[posF+3])  + invy*(FD[posF+3]+FU[posF+3]);
        d_Vk[posF+4] -= invx*(FL[posF+4]+FR[posF+4])  + invy*(FD[posF+4]+FU[posF+4]);
        
        
#if SPL>0        
        if(current_time<d_tCC[timeCC-2]){
            // Strip down
            if(j<SPL){
                // bottom left corner
                if(i<SPL){
                    spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, i, 2*npx+npy+j, i, j, d_CC, d_tCC, timeCC, current_time);
                }
                // Bottom right corner
                else if(i>npx-SPL){
                    spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, i, npx+j, npx-1-i, j, d_CC, d_tCC, timeCC, current_time);
                }
                // Central strip
                else{
                    spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, i, -1, j, -1, d_CC, d_tCC, timeCC, current_time);
                }
            }
            else{
                // Upper strip
                if(j>npy-SPL){
                    // Up left corner
                    if(i<SPL){
                        spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, npx+npy+i, 2*npx+npy+j, i, npy-1-j, d_CC, d_tCC, timeCC, current_time);
                    }
                    // Up right corner
                    else if(i>npx-SPL){
                        spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, npx+npy+i, npx+j, npx-1-i, npy-1-j, d_CC, d_tCC, timeCC, current_time);
                    }
                    // Central strip
                    else{
                        spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, npx+npy+i, -1, npy-1-j, -1, d_CC, d_tCC, timeCC, current_time);
                    }
                }
                else{
                    // left strip
                    if(i<SPL){
                        spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, 2*npx+npy+j, -1, i, -1, d_CC, d_tCC, timeCC, current_time);
                    }
                    // Right strip
                    if(i>npx-SPL){
                        spongeLayer(d_Ah, d_Ahu, d_Azb, npx, npy, d, pos, npx+j, -1, npx-1-i, -1, d_CC, d_tCC, timeCC, current_time);
                    }    
                }
            }
        }
#endif    
 



#if RK == 1
        
        coef0=0.0*(k==0) + 0.5*((k==1)) ;
        coef1=1.0-coef0;
        

        h   = coef0 * d_V0[posF]   + coef1 * d_Vk[posF];
        q.x = coef0 * d_V0[posF+1] + coef1 * d_Vk[posF+1];
        q.y = coef0 * d_V0[posF+2] + coef1 * d_Vk[posF+2];
        z   = coef0 * d_V0[posF+3] + coef1 * d_Vk[posF+3];
        hc  = coef0 * d_V0[posF+4] + coef1 * d_Vk[posF+4];


        
        state_filter(h,q,hc,epsilon_h);
        
        d_Vk[posF]   = h;
        d_Vk[posF+1] = q.x;
        d_Vk[posF+2] = q.y;
        d_Vk[posF+3] = z;
        d_Vk[posF+4] = hc;
        
        
        if(k==order-1){
            d_V0[posF]   = h;
            d_V0[posF+1] = q.x;
            d_V0[posF+2] = q.y;
            d_V0[posF+3] = z;
            d_V0[posF+4] = hc;
           
            d_V[posF]   = h;
            d_V[posF+1] = q.x;
            d_V[posF+2] = q.y;
            d_V[posF+3] = z;
            d_V[posF+4] = hc;

        }
        else{
            d_V[posF]   = h;
            d_V[posF+1] = q.x;
            d_V[posF+2] = q.y;
            d_V[posF+3] = z;
            d_V[posF+4] = hc;
        }
#else
        for (int i=0;i<5;i++){
            d_V[posF+i] = d_V0[posF+i] = d_Vk[posF+i];
        }
#endif
        
    }
}







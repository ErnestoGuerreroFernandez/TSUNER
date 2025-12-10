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



#ifndef _ARISTAS_VOLUMENES_H_
#define _ARISTAS_VOLUMENES_H_

#include <stdio.h>
#include "Edge_kernel.cu"
#include "Sim_parameters.hxx"



__device__ void interpolate_time(double *CC, double *tCC, long unsigned int timeCC, double current_time, double &h1, double &zb1, double2 &q1, int npx, int npy, int pos, double *tidesBd, float *tidesT, int tidesNPT){
    
    double dt;
    dt = tCC[1]-tCC[0];
    
    int t, tp;
    t = int((current_time-tCC[0])/dt) ;
    t*=(t>0);
    tp = t+1;
    int posCC = t*(2*npx+2*npy)*4 + pos*4;
    int posCCp = tp*(2*npx+2*npy)*4 + pos*4;
    
    double m, b;
    m = (CC[posCCp]-CC[posCC])/dt;
    b = CC[posCC] - m*tCC[t];
    h1 = m*current_time+b;
    // h1*=UNITS;


    m = (CC[posCCp+3]-CC[posCC+3])/dt;
    b = CC[posCC+3] - m*tCC[t];
    zb1 = m*current_time+b;

    h1 += zb1;

    m = (CC[posCCp+1]-CC[posCC+1])/dt;
    b = CC[posCC+1] - m*tCC[t];
    q1.x = m*current_time+b;
    q1.x *= h1   ;

    m = (CC[posCCp+2]-CC[posCC+2])/dt;
    b = CC[posCC+2] - m*tCC[t];
    q1.y = m*current_time+b;
    q1.y *= h1    ;


    // tides
    // dt = tidesT[1]-tidesT[0];
    // t = int((current_time-tidesT[0])/dt) ;
    // t*=(t>0);
    // tp = t+1;
    // posCC = t*(2*npx+2*npy) + pos;
    // posCCp = tp*(2*npx+2*npy) + pos;
    //
    // m = (tidesBd[posCCp]-tidesBd[posCC])/dt;
    // b = tidesBd[posCC] - m*tidesT[t];
    // h1 = h1 + (m*current_time+b) * (current_time<tidesT[tidesNPT-1]);
    
}


__global__ void Volumes (int *fixH, double *d_dVh, double2 *d_dVhu, double *d_dVH, double *d_Ah, double2 *d_Ahu, double *d_Azb, double *d0_dVh, double2 *d0_dVhu, double *d0_dVH, int npx, int npy, double delta_T, double dx, double dy, double current_time, double epsilon_h, double vmax, double* reconsLh, double* reconsRh, double* reconsDh, double* reconsUh, double* reconsLzb, double* reconsRzb, double* reconsDzb, double* reconsUzb, double2* reconsLhu, double2* reconsRhu, double2* reconsDhu, double2* reconsUhu, double *FL, double *FR, double *FU, double *FD, double *d_CC, double *d_tCC, long unsigned int timeCC, double *tidesBd, float *tidesT, int tidesNPT){

	int i = blockIdx.x*NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
	int j = blockIdx.y*NUM_THREADS_HORIZONTAL_VOL  + threadIdx.y;
    int pos0, pos1;
    double h0, h1, zb0, zb1;
    double2 q0, q1;
    int fixed_sediment=0; // if 0, it does not move. If 1, it moves normally.

    if ((i <= npx) && (j <= npy)) {
            if((i <= npx) && (j < npy)){
                // Riemann problems at vertical edges

                pos0 = i-1 + j*npx;
                pos1 = i + j*npx;
                
                
                if(i==0){
                    
#if WBC==2
                    
#else
                    pos0 = pos1;
                    pos1 = -1;
                        
                    h1 = h0 = reconsLh[pos0];
                    zb1 = zb0 = reconsLzb[pos0];
                    q1.x = q0.x = reconsLhu[pos0].x;
                    q1.y = q0.y = reconsLhu[pos0].y;
                        
#if WBC==1
                    q1.x = -q0.x;
#endif
#if (SPL<=0 && Read_files==1)  // Enforcing
                    if(current_time<d_tCC[timeCC-2]){
                        interpolate_time(d_CC,d_tCC,timeCC,current_time,h1,zb1,q1,npx,npy,2*npx+npy+j, tidesBd, tidesT, tidesNPT);
                    }
                    else{ // Wall boundary
                        q1.x = -q0.x;
                    }
#endif
// h1= 4.0*cos(2*PI*current_time/5.0);
// q1.x = h1*q1.x/h0;

                    fixed_sediment = fixH[pos0];
                    pos0 *=4;
                    RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FL,FR,-1,0,fixed_sediment,dx);
#endif
                }
                else{
                    if(i==npx){
#if WBC==2
                        pos1 = j*npx;
                        h0 = reconsRh[pos0];
                        zb0 = reconsRzb[pos0];
                        q0.x = reconsRhu[pos0].x;
                        q0.y = reconsRhu[pos0].y;
                        
                        h1 = reconsLh[pos1];
                        zb1 = reconsLzb[pos1];
                        q1.x = reconsLhu[pos1].x;
                        q1.y = reconsLhu[pos1].y;
                        
                        pos1 *=4;
                        
#else
                        pos1 = -2;
                        
                        h1 = h0 = reconsRh[pos0];
                        zb1 = zb0 = reconsRzb[pos0];
                        q1.x = q0.x = reconsRhu[pos0].x;
                        q1.y = q0.y = reconsRhu[pos0].y;
                        
#if WBC==1
                        q1.x = -q0.x;
#endif
                        
#endif
#if (SPL<=0 && Read_files==1)  // Enforcing
                    if(current_time<d_tCC[timeCC-2]){
                        interpolate_time(d_CC,d_tCC,timeCC,current_time,h1,zb1,q1,npx,npy,npx+j, tidesBd, tidesT, tidesNPT);
                    }
                    else{ // Wall boundary
                        q1.x = -q0.x;
                    }
#endif

                        fixed_sediment = fixH[pos0];
                        pos0 *=4;
                        RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FL,FR,1,0,fixed_sediment,dx);
                        
                    }
                    else{
                        h0 = reconsRh[pos0];
                        zb0 = reconsRzb[pos0];
                        q0.x = reconsRhu[pos0].x;
                        q0.y = reconsRhu[pos0].y;
                        
                        h1 = reconsLh[pos1];
                        zb1 = reconsLzb[pos1];
                        q1.x = reconsLhu[pos1].x;
                        q1.y = reconsLhu[pos1].y;
   
                        fixed_sediment = fixH[pos0]*fixH[pos1];
                        pos0 *=4;
                        pos1 *=4;
                        RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FL,FR,1,0,fixed_sediment,dx);
                    }
                }
                
                
                
                
            }
            if((i < npx) && (j <= npy)){
                // Riemann problems at horizontal edges
                pos0 = i + (j-1)*npx;
                pos1 = i + j*npx;
                
                
                if(j==0){
                    
#if WBC==2
                    
#else
                    pos0 = pos1;
                    pos1 = -1;
                    h1 = h0 = reconsDh[pos0];
                    zb1 = zb0 = reconsDzb[pos0];
                    q1.x = q0.x = reconsDhu[pos0].x;
                    q1.y = q0.y = reconsDhu[pos0].y;
#if WBC==1
                    q1.y = -q0.y;
#endif
                    
#if (SPL<=0 && Read_files==1)  // Enforcing
                    if(current_time<d_tCC[timeCC-2]){
                        interpolate_time(d_CC,d_tCC,timeCC,current_time,h1,zb1,q1,npx,npy,i, tidesBd, tidesT, tidesNPT);
                    }
                    else{ // Wall boundary
                        q1.y = -q0.y;
                    }
#endif

                    
                    fixed_sediment = fixH[pos0];
                    pos0 *=4;
                    RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FD,FU,0,-1,fixed_sediment,dy);
#endif
                }
                else{
                    if(j==npy){
#if WBC==2                        
                        pos1 = i;
                        
                        h0 = reconsUh[pos0];
                        zb0 = reconsUzb[pos0];
                        q0.x = reconsUhu[pos0].x;
                        q0.y = reconsUhu[pos0].y;
                        
                        h1 = reconsDh[pos1];
                        zb1 = reconsDzb[pos1];
                        q1.x = reconsDhu[pos1].x;
                        q1.y = reconsDhu[pos1].y;
                        
                        pos1*=4;
                        
#else
                        pos1 = -2;
                        h1 = h0 = reconsUh[pos0];
                        zb1 = zb0 = reconsUzb[pos0];
                        q1.x = q0.x = reconsUhu[pos0].x;
                        q1.y = q0.y = reconsUhu[pos0].y;
#if WBC==1
                        q1.y = -q0.y;
#endif
                        
#endif
#if (SPL<=0 && Read_files==1)  // Enforcing
                    if(current_time<d_tCC[timeCC-2]){
                        interpolate_time(d_CC,d_tCC,timeCC,current_time,h1,zb1,q1,npx,npy,npx+npy+i, tidesBd, tidesT, tidesNPT);
                    }
                    else{ // Wall boundary
                        q1.y = -q0.y;
                    }
#endif

// h1= 4.0*cos(2*PI*current_time/5.0);
// q1.y = h1*q1.y/h0;
                        fixed_sediment = fixH[pos0];
                        pos0 *=4;
                        RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FD,FU,0,1,fixed_sediment,dy);
                        
                    }
                    else{
                        h0 = reconsUh[pos0];
                        zb0 = reconsUzb[pos0];
                        q0.x = reconsUhu[pos0].x;
                        q0.y = reconsUhu[pos0].y;
                        
                        h1 = reconsDh[pos1];
                        zb1 = reconsDzb[pos1];
                        q1.x = reconsDhu[pos1].x;
                        q1.y = reconsDhu[pos1].y;
                        
                        
                        fixed_sediment = fixH[pos0]*fixH[pos1];
                        pos0 *=4;
                        pos1*=4;
                        RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,pos0,pos1,epsilon_h,vmax,FD,FU,0,1,fixed_sediment,dy);
                    }
                }
            }
        
    } 
    
    
    
    
}



inline void VolumeSolver(int *d_fixH, double *d_dVh, double2 *d_dVhu, double *d_dVH, double *d_Ah, double2 *d_Ahu, double *d_Azb, double *d0_dVh, double2 *d0_dVhu, double *d0_dVH, int npx, int npy, double delta_T, double dx, double dy, double current_time, double epsilon_h, double vmax, double* reconsLh, double* reconsRh, double* reconsDh, double* reconsUh, double* reconsLzb, double* reconsRzb, double* reconsDzb, double* reconsUzb, double2* reconsLhu, double2* reconsRhu, double2* reconsDhu, double2* reconsUhu, double *FL, double *FR, double *FU, double *FD, double *d_CC, double *d_tCC, long unsigned int timeCC, double *tidesBd, float *tidesT, int tidesNPT)
{
    
    
    
    dim3 blockGridEst(iDivUp(npx+1, NUM_THREADS_HORIZONTAL_VOL), iDivUp(npy+1, NUM_THREADS_VERTICAL_VOL));
    dim3 threadBlockEst(NUM_THREADS_HORIZONTAL_VOL, NUM_THREADS_VERTICAL_VOL);
    
    
    Volumes<<<blockGridEst, threadBlockEst>>>(d_fixH, d_dVh, d_dVhu, d_dVH, d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, npx, npy, delta_T, dx, dy, current_time, epsilon_h, vmax, reconsLh, reconsRh, reconsDh, reconsUh, reconsLzb, reconsRzb, reconsDzb, reconsUzb, reconsLhu, reconsRhu, reconsDhu, reconsUhu, FL, FR, FU, FD, d_CC, d_tCC, timeCC, tidesBd, tidesT, tidesNPT);

    
}
#endif 

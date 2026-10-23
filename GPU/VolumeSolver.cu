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


__global__ void VolumeSolver (int *fixH, int npx, int npy, double current_time, double epsilon_h, double* reconsL, double* reconsR, double* reconsD, double* reconsU, double *FL, double *FR, double *FU, double *FD, double *d_CC, double *d_tCC, long unsigned int timeCC, double *tidesBd, float *tidesT, int tidesNPT){

	int i = blockIdx.x*NUM_THREADS_HORIZONTAL_VOL + threadIdx.x;
	int j = blockIdx.y*NUM_THREADS_HORIZONTAL_VOL  + threadIdx.y;
    int pos0, pos1;
    int posr0, posr1;
    double h0, h1, zb0, zb1, hc0, hc1;
    double2 q0, q1;
    int fixed_sediment=0; // if 0, it does not move. If 1, it moves normally.

    if ((i <= npx) && (j <= npy)) {
            if((i <= npx) && (j < npy)){
                // Riemann problems at vertical edges

                pos0 = i-1 + j*npx;
                pos1 = i + j*npx;
                posr0 = pos0*5;
                posr1 = pos1*5;
                
                
                if(i==0){
                    
#if WBC==2
                    
#else
                    pos0 = pos1;
                    pos1 = -1;
                    posr0=posr1;
                    posr1=-1;
                        
                    h1   = h0   = reconsL[posr0];
                    q1.x = q0.x = reconsL[posr0+1];
                    q1.y = q0.y = reconsL[posr0+2];
                    zb1  = zb0  = reconsL[posr0+3];
                    hc1  = hc0  = reconsL[posr0+4];
                        
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

                    RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,posr0,posr1,epsilon_h,FL,FR,-1,0,fixed_sediment);
#endif
                }
                else{
                    if(i==npx){
#if WBC==2
                        pos1 = j*npx;
                        posr1 = pos1*5;

                        h0   = reconsR[posr0];
                        q0.x = reconsR[posr0+1];
                        q0.y = reconsR[posr0+2];
                        zb0  = reconsR[posr0+3];
                        hc0  = reconsR[posr0+4];
                        
                        h1   = reconsL[posr1];
                        q1.x = reconsL[posr1+1];
                        q1.y = reconsL[posr1+2];
                        zb1  = reconsL[posr1+3];
                        hc1  = reconsL[posr1+4];
                        

                        
#else
                        pos1 = -2;
                        posr1 = -2;
                        
                        h1   = h0   = reconsR[posr0];
                        q1.x = q0.x = reconsR[posr0+1];
                        q1.y = q0.y = reconsR[posr0+2];
                        zb1  = zb0  = reconsR[posr0+3];
                        hc1  = hc0  = reconsR[posr0+4];
                        
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
                        RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,posr0,posr1,epsilon_h,FL,FR,1,0,fixed_sediment);
                        
                    }
                    else{
                        h0   = reconsR[posr0];
                        q0.x = reconsR[posr0+1];
                        q0.y = reconsR[posr0+2];
                        zb0  = reconsR[posr0+3];
                        hc0  = reconsR[posr0+4];
                        
                        h1   = reconsL[posr1];
                        q1.x = reconsL[posr1+1];
                        q1.y = reconsL[posr1+2];
                        zb1  = reconsL[posr1+3];
                        hc1  = reconsL[posr1+4];
   
                        fixed_sediment = fixH[pos0]*fixH[pos1];
                        RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,posr0,posr1,epsilon_h,FL,FR,1,0,fixed_sediment);
                    }
                }
                
                
                
                
            }
            if((i < npx) && (j <= npy)){
                // Riemann problems at horizontal edges
                pos0 = i + (j-1)*npx;
                pos1 = i + j*npx;
                posr0 = pos0*5;
                posr1 = pos1*5;
                
                
                if(j==0){
                    
#if WBC==2
                    
#else
                    pos0 = pos1;
                    pos1 = -1;
                    posr0 = posr1;
                    posr1 = -1;

                    h1   = h0   = reconsD[posr0];
                    q1.x = q0.x = reconsD[posr0+1];
                    q1.y = q0.y = reconsD[posr0+2];
                    zb1  = zb0  = reconsD[posr0+3];
                    hc1  = hc0  = reconsD[posr0+4];
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
                    RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,posr0,posr1,epsilon_h,FD,FU,0,-1,fixed_sediment);
#endif
                }
                else{
                    if(j==npy){
#if WBC==2                        
                        pos1 = i;
                        posr1 = pos1*5;
                        
                        h0   = reconsU[posr0];
                        q0.x = reconsU[posr0+1];
                        q0.y = reconsU[posr0+2];
                        zb0  = reconsU[posr0+3]
                        hc0  = reconsU[posr0+4];
                        
                        h1   = reconsD[posr1];
                        q1.x = reconsD[posr1+1];
                        q1.y = reconsD[posr1+2];
                        zb1  = reconsD[posr1+3];
                        hc1  = reconsD[posr1+4];
                        

#else
                        pos1 = -2;
                        posr1 = -2;
                        h1   = h0   = reconsU[posr0];
                        q1.x = q0.x = reconsU[posr0+1];
                        q1.y = q0.y = reconsU[posr0+2];
                        zb1  = zb0  = reconsU[posr0+3];
                        hc1  = hc0  = reconsU[posr0+4];
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
                        RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,posr0,posr1,epsilon_h,FD,FU,0,1,fixed_sediment);
                        
                    }
                    else{
                        h0   = reconsU[posr0];
                        q0.x = reconsU[posr0+1];
                        q0.y = reconsU[posr0+2];
                        zb0  = reconsU[posr0+3];
                        hc0  = reconsU[posr0+4];
                        
                        h1   = reconsD[posr1];
                        q1.x = reconsD[posr1+1];
                        q1.y = reconsD[posr1+2];
                        zb1  = reconsD[posr1+3];
                        hc1  = reconsD[posr1+4];
                        
                        
                        fixed_sediment = fixH[pos0]*fixH[pos1];
                        RiemannSolver_Vols(h0,zb0,q0,h1,zb1,q1,posr0,posr1,epsilon_h,FD,FU,0,1,fixed_sediment);
                    }
                }
            }
        
    } 
    
    
    
    
}



#endif 

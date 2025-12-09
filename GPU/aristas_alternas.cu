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



#ifndef _ARISTAS_ALTERNAS_H_
#define _ARISTAS_ALTERNAS_H_

#include <stdio.h>
#include "Edge_kernel.cu"
#include "Sim_parameters.hxx"


void aristas_alternas(double *d_dVh, double2 *d_dVhu, double *d_dVH, double *d_Ah, double2 *d_Ahu, double *d_Azb, double *d0_dVh, double2 *d0_dVhu, double *d0_dVH, int npx, int npy, double delta_T, double dx, double dy, double current_time, double epsilon_h, double vmax, double *d_hx, double2 *d_hux, double *d_zbx, double *d_hy, double2 *d_huy, double *d_zby, int num_aristas_hor1, int num_aristas_hor2, int num_aristas_ver1, int num_aristas_ver2)
{
    
    
    // Tamano del grid y de bloque en el procesamiento de aristas
    dim3 blockGridVer1(iDivUp(num_aristas_ver1/npy, NUM_THREADS_HORIZONTAL_EDGES), iDivUp(npy, NUM_THREADS_VERTICAL_EDGES));
    dim3 blockGridVer2(iDivUp(num_aristas_ver2/npy, NUM_THREADS_HORIZONTAL_EDGES), iDivUp(npy, NUM_THREADS_VERTICAL_EDGES));
    dim3 blockGridHor1(iDivUp(npx, NUM_THREADS_HORIZONTAL_EDGES), iDivUp(num_aristas_hor1/npx, NUM_THREADS_VERTICAL_EDGES));
    dim3 blockGridHor2(iDivUp(npx, NUM_THREADS_HORIZONTAL_EDGES), iDivUp(num_aristas_hor2/npx, NUM_THREADS_VERTICAL_EDGES));

    dim3 threadBlock(NUM_THREADS_HORIZONTAL_EDGES, NUM_THREADS_VERTICAL_EDGES);
    // Tamano del grid y bloques en el calculo del nuevo estado
    // 	dim3 blockGridEst(iDivUp(npx, NUM_THREADS_HORIZONTAL_VOL), iDivUp(npy, NUM_THREADS_VERTICAL_VOL));
    // 	dim3 threadBlockEst(NUM_THREADS_HORIZONTAL_VOL, NUM_THREADS_VERTICAL_VOL);
    
    // Procesamos las aristas verticales
    procesarAristasVerGPU<<<blockGridVer1, threadBlock>>>(d_dVH, d_dVh, d_dVhu, d_hx, d_hux, d_zbx,npx, npy, delta_T, d_Ah, d_Ahu, d_Azb, dx, current_time, 1, epsilon_h, vmax, current_time);
    cudaDeviceSynchronize();
    

    procesarAristasVerGPU<<<blockGridVer2, threadBlock>>>(d_dVH, d_dVh, d_dVhu, d_hx, d_hux, d_zbx,npx, npy, delta_T, d_Ah, d_Ahu, d_Azb, dx, current_time, 2, epsilon_h,vmax, current_time);
    cudaDeviceSynchronize();
    

    
    
    //             Procesamos las aristas horizontales
    procesarAristasHorGPU<<<blockGridHor1, threadBlock>>>(d_dVH, d_dVh, d_dVhu, d_hy, d_huy, d_zby,npx, npy,delta_T, d_Ah, d_Ahu, d_Azb, dy, current_time, 1, epsilon_h,vmax, current_time);

    procesarAristasHorGPU<<<blockGridHor2, threadBlock>>>(d_dVH, d_dVh, d_dVhu, d_hy, d_huy, d_zby,npx, npy,delta_T, d_Ah, d_Ahu, d_Azb, dy, current_time, 2, epsilon_h, vmax, current_time);
//        borrar2<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, d_dVh,d_dVhu,d_dVH,npx, npy, k, epsilon_h, vmax,delta_T);
    
}
#endif

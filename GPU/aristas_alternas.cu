#ifndef _ARISTAS_ALTERNAS_H_
#define _ARISTAS_ALTERNAS_H_

#include <stdio.h>
#include "Arista_kernel.cu"
#include "Constantes.hxx"


void aristas_alternas(double *d_dVh, double2 *d_dVhu, double *d_dVH, double *d_Ah, double2 *d_Ahu, double *d_Azb, double *d0_dVh, double2 *d0_dVhu, double *d0_dVH, int npx, int npy, double delta_T, double dx, double dy, double tiempo_act, double epsilon_h, double vmax, double *d_hx, double2 *d_hux, double *d_zbx, double *d_hy, double2 *d_huy, double *d_zby, int num_aristas_hor1, int num_aristas_hor2, int num_aristas_ver1, int num_aristas_ver2)
{
    
    
    
    // Tamano del grid y de bloque en el procesamiento de aristas
    dim3 blockGridVer1(iDivUp(num_aristas_ver1/npy, NUM_HEBRAS_ANCHO_ARI), iDivUp(npy, NUM_HEBRAS_ALTO_ARI));
    dim3 blockGridVer2(iDivUp(num_aristas_ver2/npy, NUM_HEBRAS_ANCHO_ARI), iDivUp(npy, NUM_HEBRAS_ALTO_ARI));
    dim3 blockGridHor1(iDivUp(npx, NUM_HEBRAS_ANCHO_ARI), iDivUp(num_aristas_hor1/npx, NUM_HEBRAS_ALTO_ARI));
    dim3 blockGridHor2(iDivUp(npx, NUM_HEBRAS_ANCHO_ARI), iDivUp(num_aristas_hor2/npx, NUM_HEBRAS_ALTO_ARI));
    dim3 threadBlock(NUM_HEBRAS_ANCHO_ARI, NUM_HEBRAS_ALTO_ARI);
    // Tamano del grid y bloques en el calculo del nuevo estado
    // 	dim3 blockGridEst(iDivUp(npx, NUM_HEBRAS_ANCHO_EST), iDivUp(npy, NUM_HEBRAS_ALTO_EST));
    // 	dim3 threadBlockEst(NUM_HEBRAS_ANCHO_EST, NUM_HEBRAS_ALTO_EST);
    
    // Procesamos las aristas verticales
    procesarAristasVerGPU<<<blockGridVer1, threadBlock>>>(d_dVH, d_dVh, d_dVhu, d_hx, d_hux, d_zbx,npx, npy, delta_T, d_Ah, d_Ahu, d_Azb, dx, tiempo_act, 1, epsilon_h, vmax, tiempo_act);
    cudaDeviceSynchronize();
    
    
    
    
    // fprintf(stdout,"ANTES ---------------------------------------  \n");
    // borrar2<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, d_dVh,d_dVhu,d_dVH,npx, npy, k, epsilon_h, vmax,delta_T);   
    // cudaDeviceSynchronize();
    procesarAristasVerGPU<<<blockGridVer2, threadBlock>>>(d_dVH, d_dVh, d_dVhu, d_hx, d_hux, d_zbx,npx, npy, delta_T, d_Ah, d_Ahu, d_Azb, dx, tiempo_act, 2, epsilon_h,vmax, tiempo_act);
    cudaDeviceSynchronize();
    
    // fprintf(stdout,"DESPUES --------------------------------------- \n");
    // borrar2<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, d_dVh,d_dVhu,d_dVH,npx, npy, k, epsilon_h, vmax,delta_T);   
    // cudaDeviceSynchronize();
    
    
    //             Procesamos las aristas horizontales
    procesarAristasHorGPU<<<blockGridHor1, threadBlock>>>(d_dVH, d_dVh, d_dVhu, d_hy, d_huy, d_zby,npx, npy,delta_T, d_Ah, d_Ahu, d_Azb, dy, tiempo_act, 1, epsilon_h,vmax, tiempo_act);
    //            cudaDeviceSynchronize();
    procesarAristasHorGPU<<<blockGridHor2, threadBlock>>>(d_dVH, d_dVh, d_dVhu, d_hy, d_huy, d_zby,npx, npy,delta_T, d_Ah, d_Ahu, d_Azb, dy, tiempo_act, 2, epsilon_h, vmax, tiempo_act);
    //                       cudaDeviceSynchronize(); 
    
    //             aplica_fricciones<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d0_dVh, d0_dVhu,npx, npy,epsilon_h, delta_T,tiempo_act);
    //                             borrar2<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, d_dVh,d_dVhu,d_dVH,npx, npy, k, epsilon_h, vmax,delta_T);            
    //                       if(niters>1){
    //                           tiempo_act=10000;
    // 
    //                     }
    //                     niters++;
    //             borrar<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, d_dVh,d_dVhu,d_dVH,npx, npy, k, epsilon_h, vmax,delta_T);
    //             borrar2<<<blockGridEst, threadBlockEst>>>(d_Ah, d_Ahu, d_Azb, d0_dVh, d0_dVhu, d0_dVH, d_dVh,d_dVhu,d_dVH,npx, npy, k, epsilon_h, vmax,delta_T);
    
}
#endif

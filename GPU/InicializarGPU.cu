#ifndef _INICIALIZAR_GPU_H_
#define _INICIALIZAR_GPU_H_

// Función que inicializa la GPU para evitar el overhead inicial
// (ejecuta un kernel que no hace nada)
__global__ void ejecutarKernelNuloGPU() {}

extern "C" void inicializarGPU() {
	dim3 blockGrid(2, 2);
	dim3 threadBlock(16, 16);
	ejecutarKernelNuloGPU<<<blockGrid, threadBlock>>>();
}

#endif

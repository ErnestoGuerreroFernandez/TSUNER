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

/*
    This version adds multiple elements per thread sequentially.  This reduces the overall
    cost of the algorithm while keeping the work complexity O(n) and the step complexity O(log n).
    (Brent's Theorem optimization)
*/
template <class T, unsigned int blockSize>
__global__ void reduce6_min(T *g_idata, T *g_odata, unsigned int n)
{
    SharedMemory<T> smem;
    T *sdata = smem.getPointer();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;

    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridSize).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
	T thMin, val;

	if (i + blockSize >= n)
		thMin = g_idata[i];
	else
		thMin = fminf(g_idata[i], g_idata[i + blockSize]);
	i += gridSize;
	while (i < n) {
		if (i + blockSize >= n)
			val = g_idata[i];
		else
			val = fminf(g_idata[i], g_idata[i + blockSize]);
		thMin = fminf(thMin, val);
		i += gridSize;
	}
	sdata[tid] = thMin;

	__syncthreads();

    // do reduction in shared mem
    if (blockSize >= 512) {
		if (tid < 256) { sdata[tid] = thMin = fminf(sdata[tid], sdata[tid + 256]); }
		__syncthreads();
	}
    if (blockSize >= 256) {
		if (tid < 128) { sdata[tid] = thMin = fminf(sdata[tid], sdata[tid + 128]); }
		__syncthreads();
	}
    if (blockSize >= 128) {
		if (tid <  64) { sdata[tid] = thMin = fminf(sdata[tid], sdata[tid +  64]); }
		__syncthreads();
	}

/*#if (__CUDA_ARCH__ >= 300)
	if (tid < 32) {
		int warpSize = 32;

		// Fetch final intermediate sum from 2nd warp
		if (blockSize >=  64) thMin = fminf(thMin, sdata[tid + 32]);
		// Reduce final warp using shuffle
		for (int offset = warpSize/2; offset > 0; offset /= 2) {
			thMin = fminf(thMin, __shfl_down(thMin, offset));
		}
	}
#else*/
	// fully unroll reduction within a single warp
	if (tid < 32) {
		volatile T* sdatav = sdata;
		if (blockSize >=  64) sdatav[tid] = thMin = fminf(thMin, sdatav[tid + 32]);
		if (blockSize >=  32) sdatav[tid] = thMin = fminf(thMin, sdatav[tid + 16]);
		if (blockSize >=  16) sdatav[tid] = thMin = fminf(thMin, sdatav[tid +  8]);
		if (blockSize >=   8) sdatav[tid] = thMin = fminf(thMin, sdatav[tid +  4]);
		if (blockSize >=   4) sdatav[tid] = thMin = fminf(thMin, sdatav[tid +  2]);
		if (blockSize >=   2) sdatav[tid] = thMin = fminf(thMin, sdatav[tid +  1]);
	}
//#endif

	// write result for this block to global mem 
	if (tid == 0) g_odata[blockIdx.x] = thMin;
}

////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template <class T>
void reduce_min(int blocks, int threads, T *d_idata, T *d_odata, int size)
{
	dim3 dimBlock(threads, 1, 1);
	dim3 dimGrid(blocks, 1, 1);
	int smemSize = threads * sizeof(T);

	switch (threads) {
	case 512:
		reduce6_min<T, 512><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case 256:
		reduce6_min<T, 256><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case 128:
		reduce6_min<T, 128><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case 64:
		reduce6_min<T,  64><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case 32:
		reduce6_min<T,  32><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case 16:
		reduce6_min<T,  16><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case  8:
		reduce6_min<T,   8><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case  4:
		reduce6_min<T,   4><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case  2:
		reduce6_min<T,   2><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	case  1:
		reduce6_min<T,   1><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size); break;
	}
}

template
void reduce_min<float>(int blocks, int threads, float *d_idata, float *d_odata, int size);

template
void reduce_min<double>(int blocks, int threads, double *d_idata, double *d_odata, int size);

// Función que pone en blocks y threads el número de bloques y de hebras por bloque,
// respectivamente, para el kernel de reducción. Tenemos en cuenta el máximo número de
// bloques porque cada hebra puede procesar un número variable de elementos.
void getNumBlocksAndThreadsMin(int n, int maxBlocks, int maxThreads, int &blocks, int &threads)
{
/*	if (n == 1) 
		threads = 1;
	else*/
		threads = (n < maxThreads*2) ? n / 2 : maxThreads;
	blocks = n / (threads * 2);
	blocks = min(maxBlocks, blocks);
}

template <class T>
T MinReduction(T *d_data, int size)
{
	int maxBlocks = 64, maxThreads = 128;
	int numBlocks, numThreads;
	int i;
	T h_data[4096];
	T minimo;

	if (size > 4096) {
		// Run the kernel
		getNumBlocksAndThreadsMin(size, maxBlocks, maxThreads, numBlocks, numThreads);
		reduce_min<T>(numBlocks, numThreads, d_data, d_data, size);

		// We get the local mins in each GPU block until the number of elemts is equal or less than 4096
		int s = numBlocks;
		while (s > 4096) {
			getNumBlocksAndThreadsMin(s, maxBlocks, maxThreads, numBlocks, numThreads);
			reduce_min<T>(numBlocks, numThreads, d_data, d_data, s);
			s = s / (numThreads*2);
		}

		// We copy the remaining elements in CPU and finish there
		cudaMemcpy(h_data, d_data, numBlocks*sizeof(T), cudaMemcpyDeviceToHost);
		minimo = (T) 1e30;
		for (i=0; i<numBlocks; i++) {
			if (h_data[i] < minimo)
				minimo = h_data[i];
		}
	}
	else {
		// We do it in CPU
		cudaMemcpy(h_data, d_data, size*sizeof(T), cudaMemcpyDeviceToHost);
		minimo = (T) 1e30;
		for (i=0; i<size; i++) {
			if (h_data[i] < minimo)
				minimo = h_data[i];
		}
	}

	return minimo;
}


template double MinReduction<double>(double *d_in, int n);

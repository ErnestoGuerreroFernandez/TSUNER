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



#ifndef _INICIALIZAR_GPU_H_
#define _INICIALIZAR_GPU_H_

// Function that launch the GPU to avoid initial overhead. It launchs a kernel that does nothing
__global__ void ejecutarKernelNuloGPU() {}

extern "C" void initializeGPU() {
	dim3 blockGrid(2, 2);
	dim3 threadBlock(16, 16);
	ejecutarKernelNuloGPU<<<blockGrid, threadBlock>>>();
}

#endif

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



#include <stdio.h>

// This fuction returns 0 if the GPU supports CUDA, 1 if there is no GPU and 2 if there is, but it doeesn't support CUDA.
extern "C" int checkCUDAsupport(int cudaDevice)
{
    int valor, dev, deviceCount;
    valor=0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0)
        valor = 1;
	else {
        if (cudaDevice>=deviceCount) 
            valor=3;
        else {
            for (dev=0; dev < deviceCount; ++dev) {
                cudaDeviceProp deviceProp;
                cudaGetDeviceProperties(&deviceProp, dev);
                if (deviceProp.major >= 1)
                    break;
            }
            if (dev == deviceCount)
                valor = 2;
            else
                valor = 0;
        }
    }
	return valor;
}

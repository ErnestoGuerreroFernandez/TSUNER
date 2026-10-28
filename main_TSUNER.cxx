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




#include "TSUNER.h"


int main(int argc, char *argv[])
{
    double2 *dVhu;  // hux huy
    double *dVh;  // h
    double* dVH; //batimetria
    double *CC;  // Boundary conditions
    double *tCC;  // Boundary conditions time
    int *fixH; // Fixed bathymetry
    double *dyrup; // Dynamic rupture
    float *time_dyrup; // time of dynamic rupture
    long unsigned int timeCC;
    double *tides_bd;
    float *tides_time;
    int tidesNPT;
    int CUDAsupport, err;
    int cudaDevice;
    SimParams params = {};
    string input_file;
    string input_prefix;      // Prefijo de los ficheros NetCDF
    double time_gpu;

    // Fichero de datos
    input_file = argv[1];
    if (argc==3) {
        cudaDevice=atoi(argv[2]);
    }
    else {
        cudaDevice=0;
    }
    if (! FileExist(input_file)) {
        cerr << "Error: Folder not found. '" << input_file << "'" << endl;
        return EXIT_FAILURE;
    }

    // We check if the GPU exist and supports CUDA
    CUDAsupport = checkCUDAsupport(cudaDevice);
    switch(CUDAsupport) {
        case 1:
            cerr << "Error: There is no GPU" << endl;
            return EXIT_FAILURE;
            break;
        case 2:
            cerr << "Error: There is no GPU that supports CUDA" << endl;
            return EXIT_FAILURE;
            break;
        case 3:
            cerr << "N. Device erroneo" << endl;
            return EXIT_FAILURE;
            break;
    }
    initializeGPU();

    err = LoadProblemData(input_file, input_prefix, &dVh, &dVhu, &dVH, &CC, &tCC, &fixH, &dyrup, &time_dyrup, timeCC, params, &tides_bd, &tides_time, tidesNPT);
    if (err > 0)
        return EXIT_FAILURE;
    cout << scientific;
    // Sponge layer dumping function and values at the boundary
    double *d=nullptr;
    #if SPL>0
    d = new double [SPL]();
    for (int i=0; i<SPL; i++){
        d[i] = (double(SPL-i)/double(SPL))*(double(SPL-i)/double(SPL));
        //             d[i] = 1.0;
    }
    #endif

printf("%i\n",params.nvol);

    err = TSUNER(dVh, dVhu, dVH, fixH, params, (char *) input_prefix.c_str(), d, CC, tCC, timeCC, dyrup, time_dyrup, time_gpu, cudaDevice, tides_bd, tides_time, tidesNPT);
    if (err > 0) {
        if (err == 1)
            cerr << "Error: Not enough GPU memory" << endl;
        else if (err == 2)
            cerr << "Error: Not enough CPU memory" << endl;
        freeMemoryCPU(dVh,dVhu,dVH,fixH, CC, tCC, dyrup, time_dyrup);
        return EXIT_FAILURE;
    }
    cout << endl << "Total simulation time: " << time_gpu << " sec." << endl;
    freeMemoryCPU(dVh, dVhu, dVH,fixH,CC, tCC, dyrup, time_dyrup);
    #if SPL>0
    delete [] d;
    #endif
    cout<<"Normal End."<<endl;

    return EXIT_SUCCESS;
}


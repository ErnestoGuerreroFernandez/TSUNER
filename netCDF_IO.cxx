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



#include <iostream>
#include <cmath>  
#include <netcdf.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


using namespace std;


void handle_error(int val){
    std::cout<<"There has been an error in the NetCDF. The value of the error is: "<<val<<". End of the program."<<std::endl;
    exit(-1);
}


int main(){
    int ncid, varid;
    int ncid2;
    int status;

    status = nc_open("E0_3_movie2.nc", NC_WRITE, &ncid);
    if (status != NC_NOERR) handle_error(status);



    int dimId;
    long unsigned int npx, npy, time;
    nc_inq_dimid(ncid, "x", &dimId);
    nc_inq_dimlen(ncid, dimId, &npx);
    
    nc_inq_dimid(ncid, "y", &dimId);
    nc_inq_dimlen(ncid, dimId, &npy);
    
    nc_inq_dimid(ncid, "time", &dimId);
    nc_inq_dimlen(ncid, dimId, &time);

    static size_t start[] = {0,0,0};
    static size_t count[3];
    count[0] = time;
    count[1] = npy;
    count[2] = npx;
    
    double *eta, *bathy;
    eta = new double [time*npx*npy];
    bathy = new double [time*npx*npy];

    status = nc_inq_varid(ncid, "eta", &varid);
    if (status != NC_NOERR) handle_error(status);
    status = nc_get_vara_double(ncid, varid, start, count, eta);
    if (status != NC_NOERR) handle_error(status);
    
    
    status = nc_inq_varid(ncid, "b", &varid);
    if (status != NC_NOERR) handle_error(status);
    status = nc_get_vara_double(ncid, varid, start, count, bathy);
    if (status != NC_NOERR) handle_error(status);
    

    
    for (int t=0;t<time;t++){
        for (int i=0;i<npx;i++){
            for(int j=0;j<npy;j++){
                int pos = t*npx*npy + j*npx + i;
                if( abs(eta[pos]-bathy[pos]) < 1e-10 ) {
                    eta[pos] = NAN;
                }

            }
        }
    }

    status = nc_open("sim_copia.nc", NC_WRITE, &ncid2);
    if (status != NC_NOERR) handle_error(status);

    status = nc_inq_varid(ncid2, "eta", &varid);
    if (status != NC_NOERR) handle_error(status);
    
    status = nc_put_vara_double(ncid2, varid, start, count, eta);
    if (status != NC_NOERR) handle_error(status);
    
    
    status = nc_close(ncid);
    status = nc_close(ncid2);    
    delete [] eta;
    delete [] bathy;
    
}





#include <iostream>
#include <cmath>  
#include <netcdf.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// #include<detector_header.h>

using namespace std;


void handle_error(int val){
    std::cout<<"Ha habido algun error con el netCDF. El valor del error es: "<<val<<". Fin del programa."<<std::endl;
    exit(-1);
}


int main(){
    
    
    int ncid, varid;
    int ncid2;

    //Variable de error 
    int status;
    
//     strcat(file_name,"_limited.nc");

    // Abrimos archivo para edicion con NC_WRITE
    status = nc_open("sim.nc", NC_WRITE, &ncid);
    if (status != NC_NOERR) handle_error(status);
    

    // Cogemos dimensiones
    int dimId;
    long unsigned int npx, npy, time;
//     long unsigned int npx;
    nc_inq_dimid(ncid, "x", &dimId);
    nc_inq_dimlen(ncid, dimId, &npx);
    
    nc_inq_dimid(ncid, "y", &dimId);
    nc_inq_dimlen(ncid, dimId, &npy);
    
//     long unsigned int time;
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

    
    
        // Cogemos la id de la variable
    status = nc_inq_varid(ncid2, "eta", &varid);
    if (status != NC_NOERR) handle_error(status);
    
    // Y la escribimos
    status = nc_put_vara_double(ncid2, varid, start, count, eta);
    if (status != NC_NOERR) handle_error(status);
    
    
    // Cierre de ficheros
    status = nc_close(ncid);
    status = nc_close(ncid2);    
    delete [] eta;
    delete [] bathy;
    
}





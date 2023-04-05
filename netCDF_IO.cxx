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

/*

void size_netCDF(char *path, int &size, int nFiles){

    size = 0;
    char file_name[100];
    
    for (int n=0; n<nFiles; n++){
        // Variables para el netcdf
        int ncid;

        //Variable de error 
        int status;
        
        // We build netCDF file name
        strcpy(file_name,path);
        strcat(file_name,std::to_string(n).c_str());
        strcat(file_name,".nc");

        // Abrimos archivo para edicion con NC_WRITE
        status = nc_open(file_name, NC_WRITE, &ncid);
        if (status != NC_NOERR) handle_error(status);

        // Cogemos dimensiones
        int dimId;
        long unsigned int npx, t;
        nc_inq_dimid(ncid, "x", &dimId);
        nc_inq_dimlen(ncid, dimId, &npx);
        
        nc_inq_dimid(ncid, "time", &dimId);
        nc_inq_dimlen(ncid, dimId, &t);
        

        // Cierre de ficheros
        status = nc_close(ncid);
        if (status != NC_NOERR) handle_error(status);

        size += t * int(npx/NDOF);
//         *npx = npx/NDOF;
//         *time = t;
    }
}


void leer_netcdf(char *file_name, long unsigned int &npx, long unsigned int &time, double *dx, double *u, double *uL, double *uR, double *beta, double *wGP, int flag){
    
    // Variables para el netcdf
    int ncid, varid;

    //Variable de error 
    int status;
    
//     strcat(file_name,"_limited.nc");

    // Abrimos archivo para edicion con NC_WRITE
    status = nc_open(file_name, NC_WRITE, &ncid);
    if (status != NC_NOERR) handle_error(status);

    // Cogemos dimensiones
    int dimId;
//     long unsigned int npx;
    nc_inq_dimid(ncid, "x", &dimId);
    nc_inq_dimlen(ncid, dimId, &npx);
    
//     long unsigned int time;
    nc_inq_dimid(ncid, "time", &dimId);
    nc_inq_dimlen(ncid, dimId, &time);
    
//     double xL;
//     status = nc_get_att_double(ncid, NC_GLOBAL, "xL", &xL);
//     if (status != NC_NOERR) handle_error(status);
    
//     double dx;
    status = nc_get_att_double(ncid, NC_GLOBAL, "delta_x", &dx[0]);
    if (status != NC_NOERR) handle_error(status);
    
    
    // Para leer variables dependientes del tiempo de forma segura. El tamaño de la lista es tantas como dimensiones tenga el problema. Empieza en start=0 y el numero de variables que se leen son count.
    static size_t start[] = {0,0};
    static size_t count[2];
    count[0] = time;
    count[1] = npx;

    double *uaux;
    uaux    = new double[time*npx];

    // Get the varid of the data variable, based on its name.
    if (flag==0)
        status = nc_inq_varid(ncid, "v0", &varid);
    else
        status = nc_inq_varid(ncid, "v1", &varid);
    if (status != NC_NOERR) handle_error(status);
    status = nc_get_vara_double(ncid, varid, start, count, uaux);
    if (status != NC_NOERR) handle_error(status);
    
    npx = int(npx/NDOF);
    count[1] = npx;

    
    status = nc_inq_varid(ncid, "beta", &varid);
    if (status != NC_NOERR) handle_error(status);
    status = nc_get_vara_double(ncid, varid, start, count, beta);
    if (status != NC_NOERR) handle_error(status);    

    

    
    
    for (int t=0; t<time; t++){
        for (int i=0; i<npx; i++){
            int posI = t*npx + i;
            u[posI] = 0.0;
            for (int k=0; k<NDOF; k++){
                int pos = posI*NDOF + k;
                u[posI] +=  uaux[pos]*wGP[k];
            }
        }
    }
    
    // The solution is scaled in [-1,1]
    double max;
    double min;
    for (int t=0; t<time; t++){
        max = -1e30;
        min = 1e30;
        for (int i=0; i<npx; i++){
                int pos = i + t*npx ;
                if (u[pos]>max) max = u[pos];
                if (u[pos]<min) min = u[pos];
        }
        if(max==min) max += 1e-16;
        for (int i=0; i<npx; i++){
            int pos = t*npx + i;
            u[pos] = (u[pos]-min) * 2.0/(max-min) - 1.0;
        }
    }
    //
    
    
    for (int t=0; t<time; t++){
        for (int i=0; i<npx; i++){
            int pos = t*npx + i;
            int posL = pos - 1*(i>0);
            int posR = pos + 1*(i<npx-1);
            
            uL[pos] = u[posL];
            uR[pos] = u[posR];
            dx[pos] = dx[0];
        }
    }
    
    
    delete [] uaux;
    
    // Cierre de ficheros
    status = nc_close(ncid);
    if (status != NC_NOERR) handle_error(status);
    
}




void write_netcdf(char *file_name, double *p, int flag){
    
    // Variables para el netcdf
    int ncid, varid;

    //Variable de error 
    int status;
    
    int size = strlen(file_name); //Total size of string
    file_name[size-3] = '\0';  // We remove the ".nc"
    
    strcat(file_name,"_limited.nc");


    // Abrimos archivo para edicion con NC_WRITE
    status = nc_open(file_name, NC_WRITE, &ncid);
    if (status != NC_NOERR) handle_error(status);
    
    
    // Cogemos dimensiones
    int dimId;
    long unsigned int npx;
    nc_inq_dimid(ncid, "x", &dimId);
    nc_inq_dimlen(ncid, dimId, &npx);
    
    long unsigned int time;
    nc_inq_dimid(ncid, "time", &dimId);
    nc_inq_dimlen(ncid, dimId, &time);
    
    
    // Para leer variables dependientes del tiempo de forma segura. El tamaño de la lista es tantas como dimensiones tenga el problema. Empieza en start=0 y el numero de variables que se leen son count.
    static size_t start[] = {0,0};
    static size_t count[2];
    count[0] = time;
    count[1] = npx;

    
    // Cogemos la id de la variable
    if(flag==0){
        status = nc_inq_varid(ncid, "v0", &varid);
    }
    else if (flag==1){
        status = nc_inq_varid(ncid, "v1", &varid);
    }
    if (status != NC_NOERR) handle_error(status);
    
    // Y la escribimos
    status = nc_put_vara_double(ncid, varid, start, count, p);
    if (status != NC_NOERR) handle_error(status);

    
    // Cierre de ficheros
    status = nc_close(ncid);
    if (status != NC_NOERR) handle_error(status);
    

}
*/

int main(){
    
    
    int ncid, varid;
    int ncid2;

    //Variable de error 
    int status;
    
//     strcat(file_name,"_limited.nc");

    // Abrimos archivo para edicion con NC_WRITE
    status = nc_open("salida_lunes_sol.nc", NC_WRITE, &ncid);
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
// int pos = t*npx*npy + i*npy + j;
                
                if( abs(eta[pos]-bathy[pos]) < 1e-10 ) {
//                     printf("%f \n",abs(eta[pos]-bathy[pos]));
//                     cout<<pos<<endl;
                    eta[pos] = NAN;
                }

            }
        }
    }
    
    status = nc_open("salida_lunes_sol_copia.nc", NC_WRITE, &ncid2);
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





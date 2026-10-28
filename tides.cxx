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


int ncid_tides;
bool ErrorNetCDF_tides;



void check_errGRD_tides(int iret)
{
    if ((iret != NC_NOERR) && (! ErrorNetCDF_tides)) {
        fprintf(stderr, "%s\n", nc_strerror(iret));
        ErrorNetCDF_tides = true;
    }
}


void abrir_GRD_tides(const char *nombre_fich)
{
    int iret;

    iret = nc_open(nombre_fich, NC_NOWRITE, &ncid_tides);
    check_errGRD_tides(iret);
}

void cerrarGRD_tides(int id)
{
    int iret;

    iret = nc_close(id);
    check_errGRD_tides(iret);
}

void leerTamanoMallaGRD_tides(const char *nombre_fich,int *nvx, int *nvy) {
    int id, iret;
    int lon_dim, lat_dim;
    size_t nx, ny;

    iret = nc_open(nombre_fich, NC_NOWRITE, &id);
    check_errGRD_tides(iret);
    iret = nc_inq_dimid(id, "lon", &lon_dim);
    if (iret != NC_NOERR) {
        iret = nc_inq_dimid(id, "x", &lon_dim);
    }
    if (iret != NC_NOERR) {
        iret = nc_inq_dimid(id, "LON", &lon_dim);
        check_errGRD_tides(iret);
    }
    iret = nc_inq_dimid(id, "lat", &lat_dim);
    if (iret != NC_NOERR) {
        iret = nc_inq_dimid(id, "y", &lat_dim);
    }
    if (iret != NC_NOERR) {
        iret = nc_inq_dimid(id, "LAT", &lat_dim);
        check_errGRD_tides(iret);
    }
    iret = nc_inq_dimlen(id, lon_dim, &nx);
    check_errGRD_tides(iret);
    iret = nc_inq_dimlen(id, lat_dim, &ny);
    check_errGRD_tides(iret);
    iret = nc_close(id);
    check_errGRD_tides(iret);
    *nvx = nx;
    *nvy = ny;

}
void leerLongitudGRD_tides(double *lon)
{
    int id, iret;
    int lon_id;
    id=ncid_tides;

    iret = nc_inq_varid(id, "lon", &lon_id);
    if (iret != NC_NOERR) {
        iret = nc_inq_varid(id, "x", &lon_id);
    }
    if (iret != NC_NOERR) {
        iret = nc_inq_varid(id, "LON", &lon_id);
        check_errGRD_tides(iret);
    }
    iret = nc_get_var_double(id, lon_id, lon);
    check_errGRD_tides(iret);
}

void leerLatitudGRD_tides(double *lat, int flag=0)
{
    int id, iret;
    int lat_id;


    id=ncid_tides;

    // El nombre de la variable puede ser lat o y
    iret = nc_inq_varid(id, "lat", &lat_id);
    if (iret != NC_NOERR) {
        iret = nc_inq_varid(id, "y", &lat_id);
    }
    if (iret != NC_NOERR) {
        iret = nc_inq_varid(id, "LAT", &lat_id);
        check_errGRD_tides(iret);
    }
    iret = nc_get_var_double(id, lat_id, lat);
    check_errGRD_tides(iret);
}

int leerNC_time_tides(float **tides_time)
{
    int id, iret;
    int z_id;
    int dimId;
    long unsigned int num_volx, num_voly;


    long unsigned int time;
    iret = nc_inq_dimid(ncid_tides, "time", &dimId);
    if (iret != NC_NOERR) {
        iret = nc_inq_dimid(ncid_tides, "TIME", &dimId);
        check_errGRD_tides(iret);
    }
    check_errGRD_tides(iret);
    iret = nc_inq_dimlen(ncid_tides, dimId, &time);
    check_errGRD_tides(iret);

    (*tides_time) = new float[time]();

    iret = nc_get_var_float(ncid_tides, dimId, *tides_time);
    check_errGRD_tides(iret);

    return time;





}


inline double l0(double x, double xi, double xip, double xim, double dx){
    return ((x-xi)*(x-xip)/(2*dx*dx));
}
inline double l1(double x, double xi, double xip, double xim, double dx){
    return ((x-xim)*(x-xip)/(-dx*dx));
}
inline double l2(double x, double xi, double xip, double xim, double dx){
    return ((x-xim)*(x-xi)/(2*dx*dx));
}







void tides_prep(double *eta_aux, double **tides_bd, float **tides_time, int npx, int npy, string file_tides, double xmin, double xmax, double ymin, double ymax, int &npt){

    printf("Reading tides file\n");
    fflush(stdout);
    abrir_GRD_tides(file_tides.c_str());
    int npx_grd, npy_grd;
    leerTamanoMallaGRD_tides((char *) file_tides.c_str(),&npx_grd, &npy_grd);
    double *xp;
    double *yp;
    xp=new double[npx_grd];
    yp=new double[npy_grd];
    leerLongitudGRD_tides(xp);
    leerLatitudGRD_tides(yp);

    int boundary[1];
    nc_get_att_int(ncid_tides,NC_GLOBAL,"boundaries",boundary);


    npt = leerNC_time_tides(tides_time);
    double *tides_file;
    tides_file = new double[npx_grd*npy_grd*npt]();

    int varid;
    int iret = nc_inq_varid(ncid_tides, "tides", &varid);
    check_errGRD_tides(iret);
    iret=nc_get_var_double(ncid_tides,varid,tides_file);
    check_errGRD_tides(iret);

    (*tides_bd) = new double[npt*(2*npx+2*npy)];


    double dx, dy;
    dx = (xmax-xmin)/npx;
    dy = (ymax-ymin)/npy;
    double dxt, dyt;
    dxt = (xmax-xmin)/(npx_grd);
    dyt = (ymax-ymin)/(npy_grd);



    for (int t=0; t<npt; t++){
        for (int i=0; i<npx; i++){
            double x;
            x = xmin + 0.5*dx + i*dx;
            int post, postm, postp;
            post = int((x-xmin-0.5*dxt)/dxt);
            postm = post-1*(post>0);
            postp = post+1*(post<npx_grd-1);
            double xt, xtp, xtm;
            xt = xp[0] + 0.5*dxt + post*dxt;
            xtm = xt-dxt;
            xtp = xt+dxt;

            post += t*npx_grd*npy_grd;
            postm += t*npx_grd*npy_grd;
            postp += t*npx_grd*npy_grd;


            // Down
            (*tides_bd)[i + (2*npx+2*npy)*t] = tides_file[postm] * l0(x,xt,xtp,xtm,dxt)  + tides_file[post] * l1(x,xt,xtp,xtm,dxt)  + tides_file[postp] * l2(x,xt,xtp,xtm,dxt);


            // Up
            (*tides_bd)[(npx+npy) + i + (2*npx+2*npy)*t] = tides_file[postm+npx_grd*npy_grd-npx_grd] * l0(x,xt,xtp,xtm,dxt)  + tides_file[post+npx_grd*npy_grd-npx_grd] * l1(x,xt,xtp,xtm,dxt)  + tides_file[postp+npx_grd*npy_grd-npx_grd] * l2(x,xt,xtp,xtm,dxt);

        }

        for (int j=0; j<npy; j++){
            double y;
            y = ymin + 0.5*dy + j*dy;
            int post, postm, postp;
            post = int((y-ymin-0.5*dyt)/dyt);
            postm = post-1*(post>0);
            postp = post+1*(post<npy_grd-1);
            double yt, ytp, ytm;
            yt = ymin + 0.5*dyt + post*dyt;
            ytm = yt-dyt;
            ytp = yt+dyt;
            post  *= npx_grd;
            postm *= npx_grd;
            postp *= npx_grd;
            post += t*npx_grd*npy_grd;
            postm += t*npx_grd*npy_grd;
            postp += t*npx_grd*npy_grd;






            // Right
            (*tides_bd)[npx + j + (2*npx+2*npy)*t] = tides_file[postm+npx_grd-1] * l0(y,yt,ytp,ytm,dyt)  + tides_file[post+npx_grd-1] * l1(y,yt,ytp,ytm,dyt)  + tides_file[postp+npx_grd-1] * l2(y,yt,ytp,ytm,dyt);

            // Left
            (*tides_bd)[(2*npx+npy) + j + (2*npx+2*npy)*t] = tides_file[postm] * l0(y,yt,ytp,ytm,dyt)  + tides_file[post] * l1(y,yt,ytp,ytm,dyt)  + tides_file[postp] * l2(y,yt,ytp,ytm,dyt);


        }
    }



    for(int i=0; i<npt*2*(npx+npy); i++){
        if(abs((*tides_bd)[i])>1e5){
            (*tides_bd)[i] = 0.0;
        }
    }

    if (boundary[0] == 0 ){ // We read the CI
        for (int i=0; i<npx_grd*npy_grd;i++){
            eta_aux[i] = tides_file[i];
        }

    }



    delete [] tides_file;
    delete [] xp;
    delete [] yp;

}









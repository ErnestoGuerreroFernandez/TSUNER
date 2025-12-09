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



#ifndef _NETCDF_H_
#define _NETCDF_H_

#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <netcdf.h>
#include "Sim_parameters.hxx"

bool ErrorEnNetCDF;
int ncid;
int timef;
int eta_id;
int vel_id[2];
int hmax_id;
int zb_id;
int bmin_id;

void check_err(int iret)
{
	if ((iret != NC_NOERR) && (! ErrorEnNetCDF)) {
		fprintf(stderr, "%s\n", nc_strerror(iret));
		ErrorEnNetCDF = true;
	}
}


void fgennc(double xmin, double xmax, double ymin, double ymax, double* xgrid, double* ygrid, char *input_prefix, int *p_ncid, int *time_id, int& eta_id, int *vel_id, int npx, int npy, double total_time, double *bati, double R, double DS, double FR, double porosity, double OC)
{
	char file_name[256];
	char cadena[256];
	// Dimensions
// 	int grid_x_dim, grid_y_dim;
	int grid_dims[2];
	int x_dim, y_dim;
    int var_dims[3];
    int time_dim;
	// Ids
	int ncid;
// 	int grid_id, grid_x_id, grid_y_id;
	int x_id, y_id;
	double val_double, fill_double;
	struct timeval tv;
	char fecha_act[24];
	int iret;

    sprintf(file_name, "%s.nc", input_prefix);
	iret = nc_create(file_name, NC_CLOBBER, p_ncid);
	check_err(iret);
	ncid = *p_ncid;

	// Define dimensions
	iret = nc_def_dim(ncid, "x", npx, &x_dim);
	check_err(iret);
	iret = nc_def_dim(ncid, "y", npy, &y_dim);
	check_err(iret);
    
	iret = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	check_err(iret);

	// Define variables
	iret = nc_def_var(ncid, "x", NC_DOUBLE, 1, &x_dim, &x_id);
	check_err(iret);
	iret = nc_def_var(ncid, "y", NC_DOUBLE, 1, &y_dim, &y_id);
	check_err(iret);
    
    grid_dims[0] = y_dim;
    grid_dims[1] = x_dim;
    
    var_dims[0] = time_dim;
    var_dims[1] = y_dim;
    var_dims[2] = x_dim;
    iret = nc_def_var(ncid, "b", NC_DOUBLE, 3, var_dims, &zb_id);
    check_err(iret);
    
    iret=nc_def_var(ncid,"hmax",NC_DOUBLE,2,grid_dims,&hmax_id);	
	check_err(iret);
	iret=nc_def_var(ncid,"bmin",NC_DOUBLE,2,grid_dims,&bmin_id);
	check_err(iret);
    var_dims[0] = time_dim;
    var_dims[1] = y_dim;
    var_dims[2] = x_dim;
    
	grid_dims[0] = y_dim;
	grid_dims[1] = x_dim;
	iret = nc_def_var(ncid, "time", NC_DOUBLE, 1, &time_dim, time_id);
	check_err(iret);
	
    iret = nc_def_var(ncid, "eta", NC_DOUBLE, 3, var_dims, &eta_id);
    check_err(iret);
    iret = nc_put_att_text(ncid, eta_id, "units", 6, "meters");
    check_err(iret);
    iret = nc_put_att_text(ncid, eta_id, "long_name", 5, "Depth");
    check_err(iret);
    
    iret = nc_def_var(ncid, "ux", NC_DOUBLE, 3, var_dims, &(vel_id[0]));
    check_err(iret);
    iret = nc_put_att_text(ncid, vel_id[0], "units", 13, "meters/second");
    check_err(iret);
    iret = nc_put_att_text(ncid, vel_id[0], "long_name", 10, "x-velocity");
    check_err(iret);

    iret = nc_def_var(ncid, "uy", NC_DOUBLE, 3, var_dims, &(vel_id[1]));
    check_err(iret);
    iret = nc_put_att_text(ncid, vel_id[1], "units", 13, "meters/second");
    check_err(iret);
    iret = nc_put_att_text(ncid, vel_id[1], "long_name", 10, "y-velocity");
    check_err(iret);
    
            
   
	fill_double = -1e+30;
    

    iret = nc_put_att_text(ncid, hmax_id, "long_name", 4, "hmax");
    check_err(iret);
    iret = nc_put_att_text(ncid, hmax_id, "standard_name", 4, "hmax");
    check_err(iret);
    iret = nc_put_att_text(ncid, hmax_id, "units", 6, "meters");
    check_err(iret);
    iret = nc_put_att_double(ncid, hmax_id, "missing_value", NC_DOUBLE, 1, &fill_double);
    check_err(iret);
    iret = nc_put_att_double(ncid, hmax_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
    check_err(iret);

	iret = nc_put_att_text(ncid, bmin_id, "long_name", 4, "bmin");
	check_err(iret);
	iret = nc_put_att_text(ncid, bmin_id, "standard_name", 4, "bmin");
	check_err(iret);
	iret = nc_put_att_text(ncid, bmin_id, "units", 6, "meters");
	check_err(iret);
	iret = nc_put_att_double(ncid, bmin_id, "missing_value", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid, bmin_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
	check_err(iret);



	iret = nc_put_att_text(ncid, x_id, "long_name", 5, "x dim");
	check_err(iret);
	iret = nc_put_att_text(ncid, x_id, "units", 6, "meters");
	check_err(iret);

	iret = nc_put_att_text(ncid, y_id, "long_name", 5, "y dim");
	check_err(iret);
	iret = nc_put_att_text(ncid, y_id, "units", 6, "meters");
	check_err(iret);
    
	iret = nc_put_att_text(ncid, *time_id, "long_name", 4, "Time");
	check_err(iret);
	iret = nc_put_att_text(ncid, *time_id, "units", 7, "seconds");
	check_err(iret);
	iret = nc_put_att_double(ncid, eta_id, "missing_value", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid, eta_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
    
    iret = nc_put_att_double(ncid, *vel_id, "missing_value", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
	iret = nc_put_att_double(ncid, *vel_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
	check_err(iret);
    

	// Global attributes
	iret = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF-1.0");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "title", 29, "TSUNER NOAA PMEL model output");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "creator_name", 26, "Ernesto Guerrero Fernandez");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "institution", 12, "NOAA/NCTR/UW");
	check_err(iret);
	sprintf(cadena, "https://repository.library.noaa.gov/view/noaa/67438/noaa_67438_DS1.pdf");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "Paper", strlen(cadena), cadena);
	check_err(iret);
	sprintf(cadena, "https://github.com/ErnestoGuerreroFernandez");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "Github", strlen(cadena), cadena);
	check_err(iret);

	gettimeofday(&tv, NULL);
	strftime(fecha_act, 24, "%Y-%m-%d %H:%M:%S", localtime(&(tv.tv_sec)));
	iret = nc_put_att_text(ncid, NC_GLOBAL, "history", strlen(fecha_act), fecha_act);
	check_err(iret);


	val_double = total_time;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "simulation_time", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_double = R;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "Sediment_relative_density", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_double = DS;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "Sediment_mean_diameter", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_double = FR*8;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "Friction", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_double = porosity;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "Sediment_porosity", NC_DOUBLE, 1, &val_double);
	check_err(iret);
	val_double = OC;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "Sediment_critical_theta", NC_DOUBLE, 1, &val_double);
	check_err(iret);


	val_double = total_time;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "simulation_time", NC_DOUBLE, 1, &val_double);
	check_err(iret);
    
    
    
    
    
	// End definition mode
	iret = nc_enddef(ncid);
	check_err(iret);

	// Save lon
	iret = nc_put_var_double(ncid, x_id, xgrid);
	check_err(iret);
	// Save lat
	iret = nc_put_var_double(ncid, y_id, ygrid);
	check_err(iret);
    
	
}

int initNC(char *input_prefix, int npx,int npy, double xmin, double xmax, double ymin, double ymax, double total_time,  double *bati, double R, double DS, double FR, double porosity, double OC)
{
	double *xgrid, *ygrid;
    double dx,dy;
	int i;

	ErrorEnNetCDF = false;
    xgrid = new double [npx-SPL*2]();
    ygrid = new double [npy-SPL*2]();
	if (ygrid == NULL) {
		free(xgrid);
		return 1;
	}
    dx=(xmax-xmin)/npx;
    dy=(ymax-ymin)/npy;
    
    xmin+=dx*SPL;
    xmax-=dx*SPL;
    ymin+=dy*SPL;
    ymax-=dy*SPL;

    npx-=SPL*2;
    npy-=SPL*2;

	for (i=0; i<npx; i++)
		xgrid[i] = xmin+(i+0.5)*dx;
	for (i=0; i<npy; i++)
		ygrid[i] = ymin+(i+0.5)*dy;

	fgennc(xmin, xmax, ymin, ymax, xgrid, ygrid, input_prefix, &ncid, &timef, eta_id, vel_id,npx, npy, total_time, bati, R, DS, FR, porosity, OC);

	
    delete [] xgrid;
    delete [] ygrid;

	return 0;
}
/*

void writerecs(int npx, int npy, int ncid, int time_id, int var_id, int paso, double current_time, double *var)
{
	int iret;
	double t_act = current_time;
	const size_t num = paso;
	const size_t uno = 1;
	const size_t var_start[3] = {num, 0, 0};
	const size_t var_count[3] = {1, npy, npx};

	// Guardamos el tiempo
	iret = nc_put_vara_double(ncid, time_id, &num, &uno, &t_act);
	check_err(iret);

	// Guardamos la variable var
	iret = nc_put_vara_double(ncid, var_id, var_start, var_count, var);
	check_err(iret);

	iret = nc_sync(ncid);
	check_err(iret);
}
*/
void writeEtaNC(long unsigned int npx, long unsigned int npy, int numt, double current_time, double *eta)
{
    int iret;
    const size_t num = numt;
    const size_t var_start[3] = {num, 0, 0};
    const size_t var_count[3] = {1, npy, npx};

    // We save var variable
	iret = nc_put_vara_double(ncid, eta_id, var_start, var_count, eta);
	check_err(iret);
    
	iret = nc_sync(ncid);
	check_err(iret);
}



void write_time(long unsigned int npx, long unsigned int npy, int numt, double current_time)
{
    int iret;
	double t_act = current_time;
	const size_t num = numt;
	const size_t uno = 1;
    
	// Guardamos el tiempo
	iret = nc_put_vara_double(ncid, timef, &num, &uno, &t_act);
	check_err(iret);
}

void writevelNC(long unsigned int npx, long unsigned int npy, int numt, int componente, double *var)
{
    int iret;
	const size_t num = numt;
    const size_t var_start[3] = {num, 0, 0};
    const size_t var_count[3] = {1, npy, npx};
    
	
    iret = nc_put_vara_double(ncid, vel_id[componente], var_start, var_count, var);
    check_err(iret);
    
	iret = nc_sync(ncid);
	check_err(iret);
}

void writeHNC(long unsigned int npx, long unsigned int npy, int numt, double *var)
{
    int iret;
	const size_t num = numt;
    const size_t var_start[3] = {num, 0, 0};
    const size_t var_count[3] = {1, npy, npx};
    
	
    iret = nc_put_vara_double(ncid, zb_id, var_start, var_count, var);
    check_err(iret);
    
	iret = nc_sync(ncid);
	check_err(iret);
}

void write_hmax(double *hmax) {
  int iret;
  iret = nc_put_var_double(ncid, hmax_id, hmax);
  check_err(iret);
}
void write_bmin(double *bmin) {
	int iret;
	iret = nc_put_var_double(ncid, bmin_id, bmin);
	check_err(iret);
}
void closeNC()
{
	int iret;

	iret = nc_close(ncid);
	check_err(iret);
}

#endif

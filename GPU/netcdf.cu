#ifndef _NETCDF_H_
#define _NETCDF_H_

#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <netcdf.h>
#include "Constantes.hxx"

bool ErrorEnNetCDF;
// Ids de ficheros
int ncid;
// Ids de variables
int timef;
int eta_id;
int vel_id[2];
int hmax_id;
int zb_id;

void check_err(int iret)
{
	if ((iret != NC_NOERR) && (! ErrorEnNetCDF)) {
		fprintf(stderr, "%s\n", nc_strerror(iret));
		ErrorEnNetCDF = true;
	}
}



void fgennc(double xmin, double xmax, double ymin, double ymax, double* xgrid, double* ygrid, char *prefijo, int *p_ncid, int *time_id, int& eta_id, int *vel_id, int npx, int npy, double tiempo_tot, double *bati)
{
	char nombre_fich[256];
	char cadena[256];
	// Dimensiones
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

	// Creamos el fichero y entramos en modo definición

    sprintf(nombre_fich, "%s_sol.nc", prefijo);
	iret = nc_create(nombre_fich, NC_CLOBBER, p_ncid);
	check_err(iret);
	ncid = *p_ncid;

	// Definimos dimensiones
	iret = nc_def_dim(ncid, "x", npx, &x_dim);
	check_err(iret);
	iret = nc_def_dim(ncid, "y", npy, &y_dim);
	check_err(iret);
    
//     iret = nc_def_dim(ncid, "grid_x", npx, &grid_x_dim);
//     check_err(iret);
//     iret = nc_def_dim(ncid, "grid_y", npy, &grid_y_dim);
//     check_err(iret);
	iret = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	check_err(iret);

	// Definimos variables
	iret = nc_def_var(ncid, "x", NC_DOUBLE, 1, &x_dim, &x_id);
	check_err(iret);
	iret = nc_def_var(ncid, "y", NC_DOUBLE, 1, &y_dim, &y_id);
	check_err(iret);
    
//     iret = nc_def_var(ncid, "grid_x", NC_DOUBLE, 1, &grid_x_dim, &grid_x_id);
//     check_err(iret);
//     iret = nc_def_var(ncid, "grid_y", NC_DOUBLE, 1, &grid_y_dim, &grid_y_id);
//     check_err(iret);
//     grid_dims[0] = grid_y_dim;
//     grid_dims[1] = grid_x_dim;
//     iret = nc_def_var(ncid, "b", NC_DOUBLE, 2, grid_dims, &grid_id);
//     check_err(iret);
    
    var_dims[0] = time_dim;
    var_dims[1] = y_dim;
    var_dims[2] = x_dim;
    iret = nc_def_var(ncid, "b", NC_DOUBLE, 3, var_dims, &zb_id);
    check_err(iret);
    
    iret=nc_def_var(ncid,"hmax",NC_DOUBLE,2,grid_dims,&hmax_id);	
    var_dims[0] = time_dim;
    var_dims[1] = y_dim;
    var_dims[2] = x_dim;
    
	// Nota: reutilizamos el array grid_dims
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
    
            
   
	// Asignamos attributos
	fill_double = -1e+30;
    
	
//     iret = nc_put_att_text(ncid, grid_id, "long_name", 15, "Grid bathymetry");
//     check_err(iret);
//     iret = nc_put_att_text(ncid, grid_id, "standard_name", 5, "depth");
//     check_err(iret);
//     iret = nc_put_att_text(ncid, grid_id, "units", 6, "meters");
//     check_err(iret);
//     iret = nc_put_att_double(ncid, grid_id, "missing_value", NC_DOUBLE, 1, &fill_double);
//     check_err(iret);
//     iret = nc_put_att_double(ncid, grid_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
//     check_err(iret);

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




//     iret = nc_put_att_text(ncid, grid_x_id, "long_name", 6, "x Grid");
//     check_err(iret);
//     iret = nc_put_att_text(ncid, grid_x_id, "units", 6, "meters");
//     check_err(iret);
//     iret = nc_put_att_text(ncid, grid_y_id, "long_name", 6, "y Grid");
//     check_err(iret);
//     iret = nc_put_att_text(ncid, grid_y_id, "units", 6, "meters");
//     check_err(iret);
	

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
    

	// Atributos globales
	iret = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF-1.0");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "title", 27, "LavaFlow-HySEA model output");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "creator_name", 12, "MTM2015-70490");
	check_err(iret);
	iret = nc_put_att_text(ncid, NC_GLOBAL, "institution", 6, "UMA");
	check_err(iret);
	sprintf(cadena, "No paper");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "comments", strlen(cadena), cadena);
	check_err(iret);
	sprintf(cadena, "http://path.to.paper/paper.pdf");
	iret = nc_put_att_text(ncid, NC_GLOBAL, "references", strlen(cadena), cadena);
	check_err(iret);

	gettimeofday(&tv, NULL);
	strftime(fecha_act, 24, "%Y-%m-%d %H:%M:%S", localtime(&(tv.tv_sec)));
	iret = nc_put_att_text(ncid, NC_GLOBAL, "history", strlen(fecha_act), fecha_act);
	check_err(iret);


	val_double = tiempo_tot;
	iret = nc_put_att_double(ncid, NC_GLOBAL, "simulation_time", NC_DOUBLE, 1, &val_double);
	check_err(iret);
    
    
    
    
    
	// Fin del modo definición
	iret = nc_enddef(ncid);
	check_err(iret);

	// Guardamos lon
	iret = nc_put_var_double(ncid, x_id, xgrid);
	check_err(iret);
	// Guardamos lat
	iret = nc_put_var_double(ncid, y_id, ygrid);
	check_err(iret);
    

	
//     iret = nc_put_var_double(ncid, grid_x_id, xgrid);
//     check_err(iret);
//     iret = nc_put_var_double(ncid, grid_y_id, ygrid);
//     check_err(iret);
//     iret = nc_put_var_double(ncid, grid_id, bati);
//     check_err(iret);
	
}

int initNC(char *prefijo, int npx,int npy, double xmin, double xmax, double ymin, double ymax, double tiempo_tot,  double *bati)
{
	double *xgrid, *ygrid;
    double dx,dy;
	int i;

	ErrorEnNetCDF = false;
    xgrid = new double [npx]();
    ygrid = new double [npy]();
	if (ygrid == NULL) {
		free(xgrid);
		return 1;
	}
    dx=(xmax-xmin)/npx;
    dy=(ymax-ymin)/npy;

	for (i=0; i<npx; i++)
		xgrid[i] = xmin+(i+0.5)*dx;
	for (i=0; i<npy; i++)
		ygrid[i] = ymin+(i+0.5)*dy;

	fgennc(xmin, xmax, ymin, ymax, xgrid, ygrid, prefijo, &ncid, &timef, eta_id, vel_id,npx, npy, tiempo_tot, bati);

	
    delete [] xgrid;
    delete [] ygrid;

	return 0;
}
/*

void writerecs(int npx, int npy, int ncid, int time_id, int var_id, int paso, double tiempo_act, double *var)
{
	int iret;
	double t_act = tiempo_act;
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
void writeEtaNC(long unsigned int npx, long unsigned int npy, int numt, double tiempo_act, double *eta)
{
    int iret;
    const size_t num = numt;
    const size_t var_start[3] = {num, 0, 0};
    const size_t var_count[3] = {1, npy, npx};

    // Guardamos la variable var
	iret = nc_put_vara_double(ncid, eta_id, var_start, var_count, eta);
	check_err(iret);
    
	iret = nc_sync(ncid);
	check_err(iret);
}



void write_tiempo(long unsigned int npx, long unsigned int npy, int numt, double tiempo_act)
{
    int iret;
	double t_act = tiempo_act;
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
void closeNC()
{
	int iret;

	// Cerramos los ficheros
	iret = nc_close(ncid);
	check_err(iret);
}

#endif

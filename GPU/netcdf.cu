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
#include <sys/time.h>

bool ErrorEnNetCDF;
int ncid;
int timef;
int eta_id;
int vel_id[2];
int hmax_id;
int zb_id;
int bmin_id;
int eta_max_id;
int c_id;

static void check_err(int iret)
{
    if ((iret != NC_NOERR) && (!ErrorEnNetCDF)) {
        fprintf(stderr, "%s\n", nc_strerror(iret));
        ErrorEnNetCDF = true;
    }
}

static void fgennc(double xmin, double xmax, double ymin, double ymax, double* xgrid, double* ygrid, char *input_prefix, int *p_ncid, int *time_id, int& eta_id, int *vel_id, int npx, int npy, double total_time, double *bati, double R, double DS, double FR, double porosity, double OC, double Dsg)
{
    char file_name[256];
    char cadena[256];
    int grid_dims[2];
    int x_dim, y_dim;
    int var_dims[3];
    int time_dim;
    int local_ncid;
    int x_id, y_id;
    double val_double, fill_double;
    struct timeval tv;
    char fecha_act[24];
    int iret;

    sprintf(file_name, "%s.nc", input_prefix);
    iret = nc_create(file_name, NC_CLOBBER, p_ncid);
    check_err(iret);
    local_ncid = *p_ncid;

    // Define dimensions
    iret = nc_def_dim(local_ncid, "x", npx, &x_dim);
    check_err(iret);
    iret = nc_def_dim(local_ncid, "y", npy, &y_dim);
    check_err(iret);

    iret = nc_def_dim(local_ncid, "time", NC_UNLIMITED, &time_dim);
    check_err(iret);

    // Define coordinate variables
    iret = nc_def_var(local_ncid, "x", NC_DOUBLE, 1, &x_dim, &x_id);
    check_err(iret);
    iret = nc_def_var(local_ncid, "y", NC_DOUBLE, 1, &y_dim, &y_id);
    check_err(iret);

    grid_dims[0] = y_dim;
    grid_dims[1] = x_dim;

    var_dims[0] = time_dim;
    var_dims[1] = y_dim;
    var_dims[2] = x_dim;

    iret = nc_def_var(local_ncid, "b", NC_DOUBLE, 3, var_dims, &zb_id);
    check_err(iret);

    iret = nc_def_var(local_ncid, "hmax", NC_DOUBLE, 2, grid_dims, &hmax_id);
    check_err(iret);
    iret = nc_def_var(local_ncid, "bmin", NC_DOUBLE, 2, grid_dims, &bmin_id);
    check_err(iret);
    iret = nc_def_var(local_ncid, "eta_max", NC_DOUBLE, 2, grid_dims, &eta_max_id);
    check_err(iret);

    iret = nc_def_var(local_ncid, "time", NC_DOUBLE, 1, &time_dim, time_id);
    check_err(iret);

    iret = nc_def_var(local_ncid, "eta", NC_DOUBLE, 3, var_dims, &eta_id);
    check_err(iret);
    iret = nc_put_att_text(local_ncid, eta_id, "units", 6, "meters");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, eta_id, "long_name", 5, "Depth");
    check_err(iret);

    iret = nc_def_var(local_ncid, "ux", NC_DOUBLE, 3, var_dims, &(vel_id[0]));
    check_err(iret);
    iret = nc_put_att_text(local_ncid, vel_id[0], "units", 13, "meters/second");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, vel_id[0], "long_name", 10, "x-velocity");
    check_err(iret);

    iret = nc_def_var(local_ncid, "uy", NC_DOUBLE, 3, var_dims, &(vel_id[1]));
    check_err(iret);
    iret = nc_put_att_text(local_ncid, vel_id[1], "units", 13, "meters/second");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, vel_id[1], "long_name", 10, "y-velocity");
    check_err(iret);

    iret = nc_def_var(local_ncid, "c", NC_DOUBLE, 3, var_dims, &c_id);
    check_err(iret);
    iret = nc_put_att_text(local_ncid, c_id, "units", 12, "adimensional");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, c_id, "long_name", 33, "Volumetric Sediment Concentration");
    check_err(iret);

    fill_double = -1e+30;

    iret = nc_put_att_double(local_ncid, c_id, "missing_value", NC_DOUBLE, 1, &fill_double);
    check_err(iret);
    iret = nc_put_att_double(local_ncid, c_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
    check_err(iret);

    iret = nc_put_att_text(local_ncid, hmax_id, "long_name", 4, "hmax");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, hmax_id, "standard_name", 4, "hmax");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, hmax_id, "units", 6, "meters");
    check_err(iret);
    iret = nc_put_att_double(local_ncid, hmax_id, "missing_value", NC_DOUBLE, 1, &fill_double);
    check_err(iret);
    iret = nc_put_att_double(local_ncid, hmax_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
    check_err(iret);

    iret = nc_put_att_text(local_ncid, bmin_id, "long_name", 4, "bmin");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, bmin_id, "standard_name", 4, "bmin");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, bmin_id, "units", 6, "meters");
    check_err(iret);
    iret = nc_put_att_double(local_ncid, bmin_id, "missing_value", NC_DOUBLE, 1, &fill_double);
    check_err(iret);
    iret = nc_put_att_double(local_ncid, bmin_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
    check_err(iret);

    iret = nc_put_att_text(local_ncid, eta_max_id, "long_name", 7, "eta_max");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, eta_max_id, "standard_name", 7, "eta_max");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, eta_max_id, "units", 6, "meters");
    check_err(iret);
    iret = nc_put_att_double(local_ncid, eta_max_id, "missing_value", NC_DOUBLE, 1, &fill_double);
    check_err(iret);
    iret = nc_put_att_double(local_ncid, eta_max_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
    check_err(iret);

    iret = nc_put_att_text(local_ncid, x_id, "long_name", 5, "x dim");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, x_id, "units", 6, "meters");
    check_err(iret);

    iret = nc_put_att_text(local_ncid, y_id, "long_name", 5, "y dim");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, y_id, "units", 6, "meters");
    check_err(iret);

    iret = nc_put_att_text(local_ncid, *time_id, "long_name", 4, "Time");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, *time_id, "units", 7, "seconds");
    check_err(iret);

    iret = nc_put_att_double(local_ncid, eta_id, "missing_value", NC_DOUBLE, 1, &fill_double);
    check_err(iret);
    iret = nc_put_att_double(local_ncid, eta_id, "_FillValue", NC_DOUBLE, 1, &fill_double);
    check_err(iret);

    iret = nc_put_att_double(local_ncid, vel_id[0], "missing_value", NC_DOUBLE, 1, &fill_double);
    check_err(iret);
    iret = nc_put_att_double(local_ncid, vel_id[0], "_FillValue", NC_DOUBLE, 1, &fill_double);
    check_err(iret);
    iret = nc_put_att_double(local_ncid, vel_id[1], "missing_value", NC_DOUBLE, 1, &fill_double);
    check_err(iret);
    iret = nc_put_att_double(local_ncid, vel_id[1], "_FillValue", NC_DOUBLE, 1, &fill_double);
    check_err(iret);

    // Global attributes
    iret = nc_put_att_text(local_ncid, NC_GLOBAL, "Conventions", 6, "CF-1.0");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, NC_GLOBAL, "title", 29, "TSUNER NOAA PMEL model output");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, NC_GLOBAL, "creator_name", 26, "Ernesto Guerrero Fernandez");
    check_err(iret);
    iret = nc_put_att_text(local_ncid, NC_GLOBAL, "institution", 12, "NOAA/NCTR/UW");
    check_err(iret);
    sprintf(cadena, "https://repository.library.noaa.gov/view/noaa/67438/noaa_67438_DS1.pdf");
    iret = nc_put_att_text(local_ncid, NC_GLOBAL, "Paper", strlen(cadena), cadena);
    check_err(iret);
    sprintf(cadena, "https://github.com/ErnestoGuerreroFernandez");
    iret = nc_put_att_text(local_ncid, NC_GLOBAL, "Github", strlen(cadena), cadena);
    check_err(iret);

    gettimeofday(&tv, NULL);
    strftime(fecha_act, 24, "%Y-%m-%d %H:%M:%S", localtime(&(tv.tv_sec)));
    iret = nc_put_att_text(local_ncid, NC_GLOBAL, "history", strlen(fecha_act), fecha_act);
    check_err(iret);

    val_double = total_time;
    iret = nc_put_att_double(local_ncid, NC_GLOBAL, "simulation_time", NC_DOUBLE, 1, &val_double);
    check_err(iret);
    val_double = R;
    iret = nc_put_att_double(local_ncid, NC_GLOBAL, "Sediment_relative_density", NC_DOUBLE, 1, &val_double);
    check_err(iret);
    val_double = DS;
    iret = nc_put_att_double(local_ncid, NC_GLOBAL, "Sediment_mean_diameter", NC_DOUBLE, 1, &val_double);
    check_err(iret);
    val_double = FR * 8.0;
    iret = nc_put_att_double(local_ncid, NC_GLOBAL, "Friction", NC_DOUBLE, 1, &val_double);
    check_err(iret);
    val_double = porosity;
    iret = nc_put_att_double(local_ncid, NC_GLOBAL, "Sediment_porosity", NC_DOUBLE, 1, &val_double);
    check_err(iret);
    val_double = OC;
    iret = nc_put_att_double(local_ncid, NC_GLOBAL, "Sediment_critical_theta", NC_DOUBLE, 1, &val_double);
    check_err(iret);
    val_double = Dsg;
    iret = nc_put_att_double(local_ncid, NC_GLOBAL, "Geometric_mean_size_of_suspended_sediment", NC_DOUBLE, 1, &val_double);
    check_err(iret);

    // End definition mode
    iret = nc_enddef(local_ncid);
    check_err(iret);

    // Save grid coordinates
    iret = nc_put_var_double(local_ncid, x_id, xgrid);
    check_err(iret);
    iret = nc_put_var_double(local_ncid, y_id, ygrid);
    check_err(iret);
}

int initNC(char *input_prefix, double *bati, SimParams params)
{
    double *xgrid, *ygrid;
    double dx, dy, xmin, xmax, ymin, ymax;
    int npx, npy;
    int i;
    dx = params.dx;
    dy = params.dy;
    xmin = params.xmin;
    xmax = params.xmax;
    ymax = params.ymax;
    ymin = params.ymin;
    npx = params.npx;
    npy = params.npy;

    ErrorEnNetCDF = false;
    xgrid = new double [npx - SPL * 2]();
    ygrid = new double [npy - SPL * 2]();

    dx = (xmax - xmin) / npx;
    dy = (ymax - ymin) / npy;

    xmin += dx * SPL;
    xmax -= dx * SPL;
    ymin += dy * SPL;
    ymax -= dy * SPL;

    npx -= SPL * 2;
    npy -= SPL * 2;

    for (i = 0; i < npx; i++)
        xgrid[i] = xmin + (i + 0.5) * dx;
    for (i = 0; i < npy; i++)
        ygrid[i] = ymin + (i + 0.5) * dy;

    fgennc(xmin, xmax, ymin, ymax, xgrid, ygrid, input_prefix, &ncid, &timef, eta_id, vel_id, npx, npy, params.total_time, bati, params.R, params.DS, params.FR, params.porosity, params.OC, params.Dsg);

    delete [] xgrid;
    delete [] ygrid;

    return 0;
}

static void writeCNC(int npx, int npy, int numt, double *var)
{
    int iret;
    const size_t num = numt;
    const size_t var_start[3] = {num, 0, 0};
    const size_t var_count[3] = {1, static_cast<size_t>(npy), static_cast<size_t>(npx)};

    iret = nc_put_vara_double(ncid, c_id, var_start, var_count, var);
    check_err(iret);

    iret = nc_sync(ncid);
    check_err(iret);
}

void writeEtaNC(int npx, int npy, int numt, double current_time, double *eta)
{
    int iret;
    const size_t num = numt;
    const size_t var_start[3] = {num, 0, 0};
    const size_t var_count[3] = {1, static_cast<size_t>(npy), static_cast<size_t>(npx)};

    iret = nc_put_vara_double(ncid, eta_id, var_start, var_count, eta);
    check_err(iret);

    iret = nc_sync(ncid);
    check_err(iret);
}

void write_time(int npx, int npy, int numt, double current_time)
{
    int iret;
    double t_act = current_time;
    const size_t num = numt;
    const size_t uno = 1;

    iret = nc_put_vara_double(ncid, timef, &num, &uno, &t_act);
    check_err(iret);
}

void writevelNC(int npx, int npy, int numt, int componente, double *var)
{
    int iret;
    const size_t num = numt;
    const size_t var_start[3] = {num, 0, 0};
    const size_t var_count[3] = {1, static_cast<size_t>(npy), static_cast<size_t>(npx)};

    iret = nc_put_vara_double(ncid, vel_id[componente], var_start, var_count, var);
    check_err(iret);

    iret = nc_sync(ncid);
    check_err(iret);
}

void writeHNC(int npx, int npy, int numt, double *var)
{
    int iret;
    const size_t num = numt;
    const size_t var_start[3] = {num, 0, 0};
    const size_t var_count[3] = {1, static_cast<size_t>(npy), static_cast<size_t>(npx)};

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

void write_eta_max(double *eta_max) {
    int iret;
    iret = nc_put_var_double(ncid, eta_max_id, eta_max);
    check_err(iret);
}

void closeNC()
{
    int iret;
    iret = nc_close(ncid);
    check_err(iret);
}

void saveNetCDF(double *V, double *vec, int num, int npx, int npy, int nvol, double current_time)
{
    double sim_time = current_time;
    int pos;
    int posv = 0;

    write_time(npx, npy, num, sim_time);

    // eta = h - H (plane 0: h, plane 3: H/zb)
    for (int j = SPL; j < npy - SPL; j++) {
        for (int i = SPL; i < npx - SPL; i++) {
            pos = j * npx + i;
            vec[posv++] = V[pos] - V[3 * nvol + pos];
        }
    }
    writeEtaNC(npx - 2 * SPL, npy - 2 * SPL, num, sim_time, vec);

    // ux = qx / (h + EPSILON) (plane 1: qx)
    posv = 0;
    for (int j = SPL; j < npy - SPL; j++) {
        for (int i = SPL; i < npx - SPL; i++) {
            pos = j * npx + i;
            vec[posv++] = V[1 * nvol + pos] / (V[pos] + EPSILON);
        }
    }
    writevelNC(npx - 2 * SPL, npy - 2 * SPL, num, 0, vec);

    // uy = qy / (h + EPSILON) (plane 2: qy)
    posv = 0;
    for (int j = SPL; j < npy - SPL; j++) {
        for (int i = SPL; i < npx - SPL; i++) {
            pos = j * npx + i;
            vec[posv++] = V[2 * nvol + pos] / (V[pos] + EPSILON);
        }
    }
    writevelNC(npx - 2 * SPL, npy - 2 * SPL, num, 1, vec);

    // H = -zb (plane 3: zb)
    posv = 0;
    for (int j = SPL; j < npy - SPL; j++) {
        for (int i = SPL; i < npx - SPL; i++) {
            pos = j * npx + i;
            vec[posv++] = -V[3 * nvol + pos];
        }
    }
    writeHNC(npx - 2 * SPL, npy - 2 * SPL, num, vec);

    // c = hc / (h + EPSILON) (plane 4: hc)
    posv = 0;
    for (int j = SPL; j < npy - SPL; j++) {
        for (int i = SPL; i < npx - SPL; i++) {
            pos = j * npx + i;
            vec[posv++] = V[4 * nvol + pos] / (V[pos] + EPSILON);
        }
    }
    writeCNC(npx - 2 * SPL, npy - 2 * SPL, num, vec);
}

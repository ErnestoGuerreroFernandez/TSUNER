#include <sys/stat.h> 
#include <fstream>
#include <netcdf.h>
#include "cini.cxx"
#include "interpolations.cxx"
#include "tides.cxx"

int ncid_bati;
int ncid_CI;

bool ErrorNetCDF;

void check_errGRD(int iret)
{
        if ((iret != NC_NOERR) && (! ErrorNetCDF)) {
                fprintf(stderr, "%s\n", nc_strerror(iret));
                ErrorNetCDF = true;
        }
}


void abrirbatiGRD(const char *nombre_fich)
{
	int iret;

	iret = nc_open(nombre_fich, NC_NOWRITE, &ncid_bati);
	check_errGRD(iret);
}
void abrirCI_GRD(const char *nombre_fich)
{
	int iret;

	iret = nc_open(nombre_fich, NC_NOWRITE, &ncid_CI);
	check_errGRD(iret);
}

void cerrarGRD(int id)
{
	int iret;

	iret = nc_close(id);
	check_errGRD(iret);
}

void leerTamanoMallaGRD(const char *nombre_fich,int *nvx, int *nvy) {
    int id, iret;
	int lon_dim, lat_dim;
	size_t nx, ny;

	iret = nc_open(nombre_fich, NC_NOWRITE, &id);
	check_errGRD(iret);
	// Leemos el tamano de la malla.
	// El nombre de las dimensiones puede ser (lon,lat) o (x,y).
	iret = nc_inq_dimid(id, "lon", &lon_dim);
	if (iret != NC_NOERR) {
		iret = nc_inq_dimid(id, "x", &lon_dim);
	}
	if (iret != NC_NOERR) {
		iret = nc_inq_dimid(id, "LON", &lon_dim);
		check_errGRD(iret);
	}
	iret = nc_inq_dimid(id, "lat", &lat_dim);
	if (iret != NC_NOERR) {
		iret = nc_inq_dimid(id, "y", &lat_dim);
	}
	if (iret != NC_NOERR) {
		iret = nc_inq_dimid(id, "LAT", &lat_dim);
        check_errGRD(iret);
	}
	iret = nc_inq_dimlen(id, lon_dim, &nx);
	check_errGRD(iret);
	iret = nc_inq_dimlen(id, lat_dim, &ny);
	check_errGRD(iret);
	iret = nc_close(id);
	check_errGRD(iret);
	*nvx = nx;
	*nvy = ny;

}
void leerLongitudGRD(double *lon, int flag=0)
{
	int id, iret;
	int lon_id;
    if (flag==0)
        id = ncid_bati;
    else
        id=ncid_CI;
    
    // El nombre de la variable puede ser lon o x
    iret = nc_inq_varid(id, "lon", &lon_id);
	if (iret != NC_NOERR) {
		iret = nc_inq_varid(id, "x", &lon_id);
	}
    if (iret != NC_NOERR) {
		iret = nc_inq_varid(id, "LON", &lon_id);
		check_errGRD(iret);
	}
	iret = nc_get_var_double(id, lon_id, lon);
	check_errGRD(iret);
}

void leerLatitudGRD(double *lat, int flag=0)
{
	int id, iret;
	int lat_id;

    if (flag==0)
        id = ncid_bati;
    else
        id=ncid_CI;
    
    // El nombre de la variable puede ser lat o y
    iret = nc_inq_varid(id, "lat", &lat_id);
	if (iret != NC_NOERR) {
		iret = nc_inq_varid(id, "y", &lat_id);
	}
    if (iret != NC_NOERR) {
		iret = nc_inq_varid(id, "LAT", &lat_id);
		check_errGRD(iret);
	}
	iret = nc_get_var_double(id, lat_id, lat);
	check_errGRD(iret);
}
void leerBatimetriaGRD(long unsigned int num_volx, long unsigned int num_voly, double *bati)
{
	int id, iret;
	int z_id;
	const size_t start[2] = {0, 0};
	const size_t count[2] = {num_voly, num_volx};

	id = ncid_bati;
	iret = nc_inq_varid(id, "bathy", &z_id);
	check_errGRD(iret);
	iret = nc_get_vara_double(id, z_id, start, count, bati);
	check_errGRD(iret);
}


void leerCI_GRD(double *ci, long unsigned int &time)
{
	int id, iret;
	int z_id;
    int dimId;
    long unsigned int num_volx, num_voly;
    
    
    //     long unsigned int time;
    iret = nc_inq_dimid(ncid_CI, "time", &dimId);
    if (iret != NC_NOERR) {
        iret = nc_inq_dimid(ncid_CI, "TIME", &dimId);
		check_errGRD(iret);
    }
    check_errGRD(iret);
    iret = nc_inq_dimlen(ncid_CI, dimId, &time);
    check_errGRD(iret);
    
    iret = nc_inq_dimid(ncid_CI, "lat", &dimId);
    if (iret != NC_NOERR) {
        iret = nc_inq_dimid(ncid_CI, "LAT", &dimId);
		check_errGRD(iret);
    }

    iret = nc_inq_dimlen(ncid_CI, dimId, &num_voly);
    check_errGRD(iret);
    
    iret = nc_inq_dimid(ncid_CI, "lon", &dimId);
    if (iret != NC_NOERR) {
        iret = nc_inq_dimid(ncid_CI, "LON", &dimId);
		check_errGRD(iret);
    }
    iret = nc_inq_dimlen(ncid_CI, dimId, &num_volx);
    check_errGRD(iret);
    
	static size_t start[3] = {0, 0, 0};
	static size_t count[3] = {time, num_voly, num_volx};
    
    double *aux;
    aux=new double[time*num_volx*num_voly]();
    if (aux == NULL) {
		cerr << "Error: No hay memoria CPU suficiente" << endl;
		exit(-1);
	}

	id = ncid_CI;
	iret = nc_inq_varid(id, "ha_lc", &z_id);
    if (iret != NC_NOERR) {
        iret = nc_inq_varid(id, "HA", &z_id);
		check_errGRD(iret);
    }
//     iret = nc_def_var_fill(id,z_id,NC_FILL, -1e-20);
    check_errGRD(iret);
	iret = nc_get_vara_double(id, z_id, start, count, aux);
	check_errGRD(iret);
    
//     for (int t=0;t<time;t++)
    for (int i=0; i<num_volx; i++) 
        for(int j=0;j<num_voly; j++)
            {
                ci[j*num_volx+i] = aux[j*num_volx+i];
            }
    
    delete [] aux;

}


// Devuelve true si existe el fichero, false en otro caso
bool existeFichero(string fichero) 
{
	struct stat stFichInfo;
	bool existe;
	int intStat;

	// Obtenemos los atributos del fichero
	intStat = stat(fichero.c_str(), &stFichInfo);
	if (intStat == 0) {
		// Hemos obtenido los atributos del fichero. Por tanto, el fichero existe.
		existe = true;
	}
	else {
		// No hemos obtenido los atributos del fichero. Notar que esto puede
		// significar que no tenemos permiso para acceder al fichero. Para
		// hacer esta comprobacion, comprobar los valores de intStat.
		existe = false;
	}

	return existe;
}

void liberarMemoria(double *dVh,double2 *dVhu, double* dVH, int *fixH, double *CC, double *tCC, double *dyrup, float *time_dyrup)
{
	if (dVhu != NULL)			    delete [] (dVhu);
    if (dVh != NULL)			    delete [] (dVh);
    if (dVH != NULL)			    delete [] (dVH);
    if (fixH != NULL)               delete [] (fixH);

    if (CC != NULL)               delete [] (CC);
    if (tCC != NULL)               delete [] (tCC);
    if (dyrup != NULL)               delete [] (dyrup);
    if (time_dyrup != NULL)               delete [] (time_dyrup);
}




// Devuelve 0 si todo ha ido bien, 1 si no existe algun fichero, y 2 si no hay memoria CPU suficiente
int cargarDatosProblema(string fich_ent, string &prefijo, double **dVh, double2 **dVhu, double **dVH, double **CC, double **tCC, int **fixH, double **dyrup, float **time_dyrup, unsigned long int &time_grid, double &xmin, double& xmax, double &ymin, double& ymax, int &npx, int &npy, int &nvol, double& dx, double& dy, double &tiempo_tot, double &tiempo_guardar, double &CFL, double& Q, double& DS, double& FR, double& Opar, double &OC, double &HQB, double &vmax, double &epsilon_h, double **tides_bd, float **tides_time, int &tidesNPT)
{
	int i, j;
    int pos;
    string file_bat;
    string file_eta, file_fix_bathy;
    string file_tides;
    double SLR_correction;
 	char c;
    double R, porosity;
    
	// Leemos los datos del problema del fichero de entrada
	ifstream fich(fich_ent.c_str());
    

    std::string line;
    std::getline(fich, line); // Comentario
    std::getline(fich, line);  
    fich>> xmin;
    fich >> c; // Salto de linea
    std::getline(fich, line);
    fich>> xmax;
    fich >> c;
    std::getline(fich, line);
    fich>> ymin;
    fich >> c;
    std::getline(fich, line);
    fich>> ymax;
    fich >> c;
    std::getline(fich, line);
	fich >> tiempo_tot;
    fich >> c;
    std::getline(fich, line);
	fich >> tiempo_guardar;
    fich >> c;
    std::getline(fich, line);
    fich >> npx;
    fich >> c;
    std::getline(fich, line);
    fich >> npy;
    fich >> c;
    std::getline(fich, line);
	fich >> CFL;
    fich >> c;
    std::getline(fich, line);
	fich >> epsilon_h;
    fich >> c;
    std::getline(fich, line);
    fich >> R;
    fich >> c;
    std::getline(fich, line);
    fich >> DS;
    fich >> c;
    std::getline(fich, line);
    fich >> FR;
    fich >> c;
    std::getline(fich, line);
    fich >> porosity;
    fich >> c;
    std::getline(fich, line);
    fich >> OC;
    fich >> c;
    std::getline(fich, line);
    fich >> prefijo;
    fich >> c;
    std::getline(fich, line);
    fich >> file_bat;
    fich >> c;
    std::getline(fich, line);
    fich >> file_eta;
    fich >> c;
    std::getline(fich, line);
    fich >> file_fix_bathy;
    fich >> c;
    std::getline(fich, line);
    fich >> SLR_correction;
    fich >> c;
    std::getline(fich, line);
    fich >> file_tides;
	fich.close();
    
#if Manning ==0
    FR/=8.0; //Chezy law
#endif
    HQB=5.0*DS;
    Q=DS*sqrt(G*(R-1.0)*DS)/(1.0-porosity);
    Opar=FR/((R-1.0)*G*DS);
    
#if Read_files==1
    
   int read_fixH=0;
    if(file_fix_bathy.compare("0")){
        read_fixH=1;
    }
   int do_tides=0;
    if(file_tides.compare("0")){
        do_tides=1;
    }

    *dyrup = NULL;
    *time_dyrup = NULL;
    
    
    npx+=2*SPL;
    npy+=2*SPL;
    
    // Reservamos memoria
    double *bathy_interpolated, *eta_interpolated;
    double *fixH_interpolated;
    nvol = npx*npy;
    bathy_interpolated=new double[nvol]();
    eta_interpolated=new double[nvol]();
    fixH_interpolated=new double[nvol]();
    
    *fixH = new int [nvol]();
    *dVhu = new double2 [nvol]();
    *dVh = new double [nvol]();
    *dVH = new double [nvol]();
        
	if (*dVH == NULL) {
		cerr << "Error: No hay memoria CPU suficiente" << endl;
        if (bathy_interpolated != NULL) delete[] bathy_interpolated;
        if (eta_interpolated != NULL) delete[] eta_interpolated;
        if (fixH_interpolated != NULL) delete[] fixH_interpolated;
		return 2;
	}
	

    // ********************   Bathy 
    printf("Reading bathymetry file\n");
    fflush(stdout);
    abrirbatiGRD((char *) file_bat.c_str());
    int npx_grd, npy_grd;
    leerTamanoMallaGRD((char *) file_bat.c_str(),&npx_grd, &npy_grd);
    double *xp;
    double *yp;
    xp=new double[npx_grd];
    yp=new double[npy_grd];
    leerLongitudGRD(xp);
    leerLatitudGRD(yp);
    
    double *bathy_aux;
    bathy_aux=new double[npx_grd*npy_grd];
    if (bathy_aux == NULL) {
        cerr << "Error: No hay memoria CPU suficiente" << endl;
        delete[] bathy_aux;
        return 2;
    }
    leerBatimetriaGRD(npx_grd, npy_grd, bathy_aux);  // Bathymetry read from file, no interpolation.
    double dx_grid, dy_grid;
    dx_grid = (xp[npx_grd-1]-xp[0])/(npx_grd-1);
    dy_grid = (yp[npy_grd-1]-yp[0])/(npy_grd-1);
    xmin = xp[0] - 0.5*dx_grid;
    xmax = xp[npx_grd-1] + 0.5*dx_grid;
    ymin = yp[0] - 0.5*dy_grid;
    ymax = yp[npy_grd-1] + 0.5*dy_grid;
    double xmin_grid, ymin_grid;
    xmin_grid=xmin;
    ymin_grid=ymin;


    
    xmin-=dx_grid*SPL;
    xmax+=dx_grid*SPL;
    ymin-=dy_grid*SPL;
    ymax+=dy_grid*SPL;
    
    interpolate_spherics(xmin,xmax,ymin,ymax, xmin_grid,ymin_grid, npx_grd, npy_grd, dx_grid, dy_grid, npx, npy, bathy_aux, bathy_interpolated);
    delete[] bathy_aux;
    delete[] xp;
    delete[] yp;
    cerrarGRD(ncid_bati);
//     printf("Done\n");
    printf("\n");
    
    // ********************   Fixed Bathy  
    if(read_fixH){
        printf("Reading fixed bathymetry file.\n");
        fflush(stdout);
        abrirbatiGRD((char *) file_fix_bathy.c_str());
        leerTamanoMallaGRD((char *) file_fix_bathy.c_str(),&npx_grd, &npy_grd);
        xp=new double[npx_grd];
        yp=new double[npy_grd];
        leerLongitudGRD(xp);
        leerLatitudGRD(yp);
        
        double *fixH_aux;
        fixH_aux=new double[npx_grd*npy_grd];
        if (fixH_aux == NULL) {
            cerr << "Error: No hay memoria CPU suficiente" << endl;
            delete[] fixH_aux;
            liberarMemoria(*dVh,*dVhu,*dVH,*fixH,*CC,*tCC,*dyrup,*time_dyrup);
            return 2;
        }
        leerBatimetriaGRD(npx_grd, npy_grd, fixH_aux);  // Bathymetry read from file, no interpolation.
        dx_grid = (xp[npx_grd-1]-xp[0])/(npx_grd-1);
        dy_grid = (yp[npy_grd-1]-yp[0])/(npy_grd-1);
        
        interpolate_spherics(xmin,xmax,ymin,ymax, xp[0], yp[0], npx_grd, npy_grd, dx_grid, dy_grid, npx, npy, fixH_aux, fixH_interpolated);
        delete[] fixH_aux;
        delete[] xp;
        delete[] yp;
        cerrarGRD(ncid_bati);
//         printf("Done\n");
        printf("\n");
    }
    
    // ********************   CI 
    printf("Reading initial conditions\n");
    fflush(stdout);
    abrirCI_GRD(file_eta.c_str());
    int npx_grd_CI, npy_grd_CI;
    leerTamanoMallaGRD((char *) file_eta.c_str(),&npx_grd_CI, &npy_grd_CI);
    double *xp_CI;
    double *yp_CI;
    xp_CI=new double[npx_grd_CI];
    yp_CI=new double[npy_grd_CI];
    leerLongitudGRD(xp_CI,1);
    leerLatitudGRD(yp_CI,1);
    double *eta_aux;
    eta_aux=new double[npx_grd_CI*npy_grd_CI];
    if (eta_aux == NULL) {
        cerr << "Error: No hay memoria CPU suficiente" << endl;
        delete[] eta_aux;
        return 2;
    }
    leerCI_GRD(eta_aux, time_grid);  // Free surface read from file, no interpolation.

    if(do_tides){
        tides_prep(eta_aux,tides_bd,tides_time, npx, npy, file_tides, xmin, xmax, ymin, ymax, tidesNPT);
    }
    else{
        tidesNPT=-1;
    }

    double dx_CI, dy_CI;
    dx_CI = (xp_CI[npx_grd_CI-1]-xp_CI[0])/npx_grd_CI;
    dy_CI = (yp_CI[npy_grd_CI-1]-yp_CI[0])/npy_grd_CI;
    interpolate_spherics(xmin,xmax,ymin,ymax, xp_CI[0], yp_CI[0], npx_grd_CI, npy_grd_CI, dx_CI, dy_CI, npx, npy, eta_aux, eta_interpolated);
    fix_eta(eta_interpolated, bathy_interpolated, npx, npy, SLR_correction);
    delete [] eta_aux;
//     printf("Done\n");
    printf("\n");

    
    // ********************** CC
    printf("Reading boundary conditions\n");
    fflush(stdout);
    *tCC=new double [time_grid]();
    *CC = new double [(2*npx+2*npy)*4*time_grid]();
    if (*CC == NULL) {
        cerr << "Error: No hay memoria CPU suficiente" << endl;
        return 2;
    }



    interpolate_CC(file_eta, (*CC), (*tCC), npx, npy, xmin, xmax, ymin, ymax,  xp_CI[0], yp_CI[0],  time_grid, dx_CI, dy_CI, bathy_interpolated, ncid_CI, SLR_correction);
    delete [] xp_CI;
    delete [] yp_CI;
//     printf("Done\n");
    printf("\n");
    
    dx = (xmax-xmin)/(npx);
    dy = (ymax-ymin)/npy;
    // ******************************* Pasamos de esfericas a cartesianas.


    double dyE, *dxE, dy_unidad, dmax;
    dxE = new double[npy]();
    
    dy_transform(dy_unidad);
    dx_transform(ymin,dy,npy, xmax, xmin, dxE, dx, dmax);
    
    xmin = 0; 
    xmax = dmax;
    ymax = dy_unidad*(ymax-ymin);
    ymin = 0;

    printf("Converting into spherical coodinates\n");
    fflush(stdout);
    dyE = dy*dy_unidad;
    // ******************************* Bathy
    double *bathy_interpolated_spheric;
    bathy_interpolated_spheric = new double[nvol]();
    spheric2cartesian(xmin, xmax, ymin, ymax, dxE, dyE, npx, npy, bathy_interpolated, bathy_interpolated_spheric);
    delete [] bathy_interpolated;
    
    // ******************************* Fix bathy
    double *fixH_interpolated_spheric;
    if(read_fixH){
        fixH_interpolated_spheric = new double[nvol]();
        spheric2cartesian(xmin, xmax, ymin, ymax, dxE, dyE, npx, npy, fixH_interpolated, fixH_interpolated_spheric);
        delete [] fixH_interpolated;
    }
    // ******************************* eta
    double *eta_interpolated_spheric;
    eta_interpolated_spheric = new double[nvol]();
    spheric2cartesian(xmin, xmax, ymin, ymax, dxE, dyE, npx, npy, eta_interpolated, eta_interpolated_spheric);
    delete [] eta_interpolated;

    // ******************************* CC
    spheric2cartesian_CC(xmin, xmax, ymin, ymax, dxE, dyE, npx, npy, time_grid, (*CC));
    delete[] dxE;
    cerrarGRD(ncid_CI);
//     printf("Done\n");
    printf("\n");

    
    
    dx = (xmax-xmin)/(npx);
    dy = (ymax-ymin)/npy;
    
    printf("Value assignation\n");
    fflush(stdout);
    for (j=0; j<npy; j++) {
        printProgress(float(j+1)/npy);
        for (i=0; i<npx; i++) {
            pos=j*npx+i;
            
            (*dVH)[pos] = (double) (bathy_interpolated_spheric[pos]) ;
            (*dVh)[pos] = (double) (eta_interpolated_spheric[pos]) *UNITS  + (*dVH)[pos]  + SLR_correction;
            // (*dVh)[pos] = 5  + (*dVH)[pos]  + SLR_correction;
            (*dVh)[pos] *= ((*dVh)[pos]>=0.0);
            (*dVhu)[pos].x = 0.0;
            (*dVhu)[pos].y = 0.0;            
            
            if(read_fixH){
                if(fixH_interpolated_spheric[pos]>=0.5){
                    (*fixH)[pos] = 1;
                }
                else{
                    (*fixH)[pos] = 0;
                }
            }
            else{
                (*fixH)[pos] = 1;
            }
            
        }   
    }

    delete[] bathy_interpolated_spheric;
    delete[] eta_interpolated_spheric;
//     printf("Done\n");
    printf("\n");

#else  // Del Read_files==1
    
    dx=(xmax-xmin)/(npx);
    dy=(ymax-ymin)/(npy);
    

    npx+=2*SPL;
    npy+=2*SPL;
    xmin-=SPL*dx;
    ymin-=SPL*dy;
    xmax+=SPL*dx;
    ymax+=SPL*dy;

    nvol = npx*npy;

    vmax = VMAX;
    
    *fixH = new int[nvol]();
    *dVhu = new double2 [nvol]();
    *dVh = new double [nvol]();
    *dVH = new double [nvol]();

    *time_dyrup = new float[40]; // It assumes range(0,400,10)
    *dyrup = new double[nvol*40]();
    
	if (*dVH == NULL) {
		cerr << "Error: No hay memoria CPU suficiente" << endl;
 		liberarMemoria(*dVh,*dVhu,*dVH,*fixH,*CC,*tCC,*dyrup,*time_dyrup);
		return 2;
	}

	for (int t=0; t<40; t++) (*time_dyrup)[t] = t*10;

    for (j=0; j<npy; j++) {
        for (i=0; i<npx; i++) {
            pos=j*npx+i;
            double x,y;
            x = xmin + 0.5*dx + i*dx;
            y = ymin + 0.5*dy + j*dy;
            (*dVh)[pos] = h_inicial(x,y);
            (*dVH)[pos]=profundidad(x,y);
            velocidad(x,y,(*dVhu)[pos]);
            (*dVhu)[pos].x *= (*dVh)[pos];
            (*dVhu)[pos].y *= (*dVh)[pos];
            (*fixH)[pos] = fixed_bed_CI(x,y);

            if(i>npx/2-5 && i<npx/2+5 && j<npy/2+5 && j>npy/2-5){
                for (int t=0;t<40;t++){
                    int posdy = pos + t *nvol;
                    (*dyrup)[posdy] = -t*0.5   +  (*dVH)[pos];  // Fijate que el valor que entra es bathymetria + perturbacion
                }
            }


        }         
    }


#endif

	return 0;
}


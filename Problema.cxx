#include <sys/stat.h> 
#include <fstream>
#include <netcdf.h>
#include "cini.cxx"


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
	// Leemos el tamaño de la malla.
	// El nombre de las dimensiones puede ser (lon,lat) o (x,y).
	iret = nc_inq_dimid(id, "lon", &lon_dim);
	if (iret != NC_NOERR) {
		iret = nc_inq_dimid(id, "x", &lon_dim);
		check_errGRD(iret);
	}
	iret = nc_inq_dimid(id, "lat", &lat_dim);
	if (iret != NC_NOERR) {
		iret = nc_inq_dimid(id, "y", &lat_dim);
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
    nc_inq_dimid(ncid_CI, "time", &dimId);
    nc_inq_dimlen(ncid_CI, dimId, &time);
    
    nc_inq_dimid(ncid_CI, "lat", &dimId);
    nc_inq_dimlen(ncid_CI, dimId, &num_voly);
    
    nc_inq_dimid(ncid_CI, "lon", &dimId);
    nc_inq_dimlen(ncid_CI, dimId, &num_volx);
    
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
	check_errGRD(iret);
	iret = nc_get_vara_double(id, z_id, start, count, aux);
	check_errGRD(iret);
    
//     for (int t=0;t<time;t++)
    for (int i=0; i<num_volx; i++) 
        for(int j=0;j<num_voly; j++)
            {
                ci[j*num_volx+i] = aux[ j*num_volx+i];
//                 printf("%.15f  \n",ci[j*num_volx+i]);
            }
    
    delete [] aux;
//     exit(-1);
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
		// hacer esta comprobaci�n, comprobar los valores de intStat.
		existe = false;
	}

	return existe;
}

void liberarMemoria(double *dVh,double2 *dVhu, double* dVH)
{
	if (dVhu != NULL)			    delete [] (dVhu);
    if (dVh != NULL)			    delete [] (dVh);
    if (dVH != NULL)			    delete [] (dVH);
}


inline double L0(double x){
    return x/2.0*(x-1);    
}
inline double L1(double x){
    return 1-x*x;    
}
inline double L2(double x){
    return x/2.0*(x+1);    
}

inline double dy_transform(double phi1, double phi2){
    return 2*Re*(phi2-phi1)/2.0 *PI/180;
}

inline double dx_transform(double alpha1, double alpha2, double phi){
    return 2*Re*asin( sqrt( cos(phi*PI/180)*cos(phi*PI/180)*sin(0.5*(alpha2-alpha1)*PI/180) *sin(0.5*(alpha2-alpha1)*PI/180) )  );
}



void interpolate_fun(double *bathy, double *p, int npx, int npy, int npx_grd, int npy_grd, double xmin, double xmax, double ymin, double ymax){
    
    double x,y;
    double dx,dy;
    double dx_grid, dy_grid;
    double xt, yt;
    int i_grid, j_grid;
    double N1,N2,N3,N4,N5,N6,N7,N8,N9;
    
    dx = (xmax-xmin)/npx;
    dy = (ymax-ymin)/npy;
    

    dx_grid = (xmax-xmin)/npx_grd;
    dy_grid = (ymax-ymin)/npy_grd;
    
    for (int j=0; j<npy; j++){
        for (int i=0; i<npx; i++){
            int pos, posp;
            int posL, posR, posU, posD, posLU, posLD, posRU, posRD;
            
            x = xmin + 0.5*dx + i*dx;
            y = ymin + 0.5*dy + j*dy;
        
            xt = 1.0/(xmax-xmin)*(2*x-(xmax+xmin));
            yt = 1.0/(ymax-ymin)*(2*y-(ymax+ymin));
        

            i_grid = int((x - xmin - 0.5*dx_grid)/dx_grid);
            j_grid = int((y - ymin - 0.5*dy_grid)/dy_grid);
            
            pos = j_grid*npx_grd + i_grid;
            posp = j*npx +i;
        
            posL = pos-1;
            posR = pos+1;
            posU = pos+npx_grd;
            posD = pos-npx_grd;
            posLU = posU-1;
            posLD = posD-1;
            posRU = posU+1;
            posRD = posD+1;
        
            if(i_grid==0){
                posL=pos;
                posLU=pos;
                posLD=pos;
            }
            if(i_grid==npx_grd-1){
                posR=pos;
                posRD=pos;
                posRU=pos;
            }
                
            if(j_grid==0){
                posD=pos;
                posLD=pos;
                posRD=pos;
            }
                
            if(j_grid==npy_grd-1){
                posU=pos;
                posLU=pos;
                posRU=pos;
            }
                
        
            N1 = L0(xt)*L0(yt);
            N2 = L1(xt)*L0(yt);
            N3 = L2(xt)*L0(yt);
            N4 = L0(xt)*L1(yt);
            N5 = L1(xt)*L1(yt);
            N6 = L2(xt)*L1(yt);
            N7 = L0(xt)*L2(yt);
            N8 = L1(xt)*L2(yt);
            N9 = L2(xt)*L2(yt);        

        
            p[posp] = bathy[posLD]*N1 + bathy[posD]*N2 + bathy[posRD]*N3 + bathy[posL]*N4 + bathy[pos]*N5 + bathy[posR]*N6 + bathy[posLU]*N7 + bathy[posU]*N8 + bathy[posRU]*N9; 
            
        }    
    }
}


void interpolate_CI(double *bathy, double *p, int npx, int npy, int npx_grd, int npy_grd, double xmin, double xmax, double ymin, double ymax, double xmin_CI, double xmax_CI, double ymin_CI, double ymax_CI, int npx_grd_CI, int npy_grd_CI){
    
    double x,y;
    double dx,dy;
    double dx_grid, dy_grid;
    double xt, yt;
    int i_grid, j_grid;
    double N1,N2,N3,N4,N5,N6,N7,N8,N9;
    
    dx = (xmax-xmin)/npx;
    dy = (ymax-ymin)/npy;
    
    dx_grid = (xmax_CI-xmin_CI)/npx_grd_CI;
    dy_grid = (ymax_CI-ymin_CI)/npy_grd_CI;
    
    for (int j=0; j<npy; j++){
        for (int i=0; i<npx; i++){
            int pos, posp;
            int posL, posR, posU, posD, posLU, posLD, posRU, posRD;
            
            x = xmin + 0.5*dx + i*dx;
            y = ymin + 0.5*dy + j*dy;
        
            xt = 1.0/(xmax-xmin)*(2*x-(xmax+xmin));
            yt = 1.0/(ymax-ymin)*(2*y-(ymax+ymin));
        

            i_grid = int((x - xmin_CI - 0.5*dx_grid)/dx_grid);
            j_grid = int((y - ymin_CI - 0.5*dy_grid)/dy_grid);
            
            pos = j_grid*npx_grd + i_grid;
            posp = j*npx +i;
        
            posL = pos-1;
            posR = pos+1;
            posU = pos+npx_grd;
            posD = pos-npx_grd;
            posLU = posU-1;
            posLD = posD-1;
            posRU = posU+1;
            posRD = posD+1;
        
            if(i_grid==0){
                posL=pos;
                posLU=pos;
                posLD=pos;
            }
            if(i_grid==npx_grd-1){
                posR=pos;
                posRD=pos;
                posRU=pos;
            }
                
            if(j_grid==0){
                posD=pos;
                posLD=pos;
                posRD=pos;
            }
                
            if(j_grid==npy_grd-1){
                posU=pos;
                posLU=pos;
                posRU=pos;
            }
                
        
            N1 = L0(xt)*L0(yt);
            N2 = L1(xt)*L0(yt);
            N3 = L2(xt)*L0(yt);
            N4 = L0(xt)*L1(yt);
            N5 = L1(xt)*L1(yt);
            N6 = L2(xt)*L1(yt);
            N7 = L0(xt)*L2(yt);
            N8 = L1(xt)*L2(yt);
            N9 = L2(xt)*L2(yt);        

        
            p[posp] = bathy[posLD]*N1 + bathy[posD]*N2 + bathy[posRD]*N3 + bathy[posL]*N4 + bathy[pos]*N5 + bathy[posR]*N6 + bathy[posLU]*N7 + bathy[posU]*N8 + bathy[posRU]*N9; 
            
        }    
    }
}


void interpolate_CC(string file, double *CC,  double *tCC, int npx, int npy, double xmin, double xmax, double ymin, double ymax, long unsigned int time, double xmin_CI, double xmax_CI, double ymin_CI, double ymax_CI, double *H){
    
    
    double x,y;
    double dx,dy;
    double dx_grid, dy_grid;
    double xt, yt;
    int i_grid, j_grid;
    int idu;
    int var_id, dimId;
	int id = ncid_CI;
    long unsigned int npx_grd, npy_grd;
    

    double *fh, *fhu, *fhv;
    

	
    // We build the file XXXu.nc
    file.erase(file.end()-4);
    file.insert(file.end()-3,'u');
	
	int iret = nc_open(file.c_str(), NC_NOWRITE, &idu);

    nc_inq_dimid(idu, "lat", &dimId);
    nc_inq_dimlen(idu, dimId, &npy_grd);
    
    nc_inq_dimid(idu, "lon", &dimId);
    nc_inq_dimlen(idu, dimId, &npx_grd);
    
    nc_inq_dimid(idu, "time", &dimId);
//     nc_inq_dimlen(idu, dimId, &tCC);
    
    static size_t start_t[1] = {0};
	static size_t count_t[1] = {time};
    
    iret = nc_inq_varid(idu, "time", &var_id);
	check_errGRD(iret);
	iret = nc_get_vara_double(idu, var_id, start_t, count_t, tCC);
	check_errGRD(iret);
    
    
    static size_t start[3] = {0, 0, 0};
	static size_t count[3] = {time, npy_grd, npx_grd};
    
    fhu=new double[time*npx_grd*npy_grd]();
    fhv=new double[time*npx_grd*npy_grd]();
    fh=new double[time*npx_grd*npy_grd]();
    if (fh == NULL) {
		cerr << "Error: No hay memoria CPU suficiente" << endl;
		exit(-1);
	}
    
    
    iret = nc_inq_varid(idu, "ua_lc", &var_id);
	check_errGRD(iret);
	iret = nc_get_vara_double(idu, var_id, start, count, fhu);
	check_errGRD(iret);
    iret = nc_close(idu);
	check_errGRD(iret);
    
    // We build the file XXXv.nc
    file.erase(file.end()-4);
    file.insert(file.end()-3,'v');
	
	iret = nc_open(file.c_str(), NC_NOWRITE, &idu);
	check_errGRD(iret);
    iret = nc_inq_varid(idu, "va_lc", &var_id);
	check_errGRD(iret);
	iret = nc_get_vara_double(idu, var_id, start, count, fhv);
	check_errGRD(iret);
    iret = nc_close(idu);
	check_errGRD(iret);

    
	iret = nc_inq_varid(id, "ha_lc", &var_id);
	check_errGRD(iret);
	iret = nc_get_vara_double(id, var_id, start, count, fh);
	check_errGRD(iret);

//     for (int i=0; i<npx_grd; i++)
//     printf("%d   %.15f  \n",i,fh[i]);
//                 printf("%d   %d  %d %.15f  %.15f  %.15f     \n",i,i_grid, j_grid, y, ymin_CI, dy_grid);
    
    dx = (xmax-xmin)/npx;
    dy = (ymax-ymin)/npy;
    
    dx_grid = (xmax_CI-xmin_CI)/npx_grd;
    dy_grid = (ymax_CI-ymin_CI)/npy_grd;
    
    for (int t=0; t<time; t++){
        int pos, pos0, pos1;
        int posC = t*(2*npx+2*npy)*4;
        double N0, N1, N2;
        for (int i=0; i<npx; i++){
            x = xmin + 0.5*dx + i*dx;
            y = ymin + 0.5*dy *0;
            xt = 1.0/(xmax-xmin)*(2*x-(xmax+xmin));
            
            N0 = L0(xt);
            N1 = L1(xt);
            N2 = L2(xt);
            
            i_grid = int((x - xmin_CI - 0.5*dx_grid)/dx_grid);
            j_grid = int((y - ymin_CI - 0.5*dy_grid)/dy_grid);
            
            // Frontera inf
            pos = t*npx_grd*npy_grd + npx_grd*j_grid + i_grid;
            pos0 = pos-1*(i_grid>0);
            pos1 = pos+1*(i_grid<npx_grd-1);
            
            CC[posC+i*4] = N0 * fh[pos0] + N1*fh[pos] + N2*fh[pos1];
            CC[posC+i*4+1] = N0 * fhu[pos0] + N1*fhu[pos] + N2*fhu[pos1];
            CC[posC+i*4+2] = N0 * fhv[pos0] + N1*fhv[pos] + N2*fhv[pos1];
            CC[posC+i*4+3] = H[i];
            
            // Frontera sup
            y = ymax - 0.5*dy *0;
            j_grid = int((y - ymin_CI - 0.5*dy_grid)/dy_grid);
            pos = t*npx_grd*npy_grd + (j_grid)*npx_grd + i_grid;
            pos0 = pos-1*(i_grid>0);
            pos1 = pos+1*(i_grid<npx_grd-1);
            
            
            CC[posC+(npx+npy+i)*4]   = N0*fh[pos0] + N1*fh[pos] + N2*fh[pos1];
            CC[posC+(npx+npy+i)*4+1] = N0*fhu[pos0] + N1*fhu[pos] + N2*fhu[pos1];
            CC[posC+(npx+npy+i)*4+2] = N0*fhv[pos0] + N1*fhv[pos] + N2*fhv[pos1];
            CC[posC+(npx+npy+i)*4+3] = H[(npy-1)*npx+i];
            

        }
        for (int j=0; j<npy; j++){
            y = ymin + 0.5*dy + j*dy;
            x = xmin + 0.5*dx *0;
            yt = 1.0/(ymax-ymin)*(2*y-(ymax+ymin));
            
            i_grid = int((x - xmin_CI - 0.5*dx_grid)/dx_grid);
            j_grid = int((y - ymin_CI - 0.5*dy_grid)/dy_grid);
            
            // Frontera izq
            pos = t*npx_grd*npy_grd + j_grid*npx_grd + i_grid;
            pos0 = pos-npx_grd*(j_grid>0);
            pos1 = pos+npx_grd*(j_grid<npy_grd-1);
            
            CC[posC+(2*npx+npy+j)*4]   = N0*fh[pos0] + N1*fh[pos] + N2*fh[pos1];
            CC[posC+(2*npx+npy+j)*4+1] = N0*fhu[pos0] + N1*fhu[pos] + N2*fhu[pos1];
            CC[posC+(2*npx+npy+j)*4+2] = N0*fhv[pos0] + N1*fhv[pos] + N2*fhv[pos1];
            CC[posC+(2*npx+npy+j)*4+3] = H[j*npx];
            
            // Frontera drch
            x = xmax - 0.5*dx *0;
            i_grid = int((x - xmin_CI - 0.5*dx_grid)/dx_grid);
            pos = t*npx_grd*npy_grd + j_grid*npx_grd + i_grid;
            pos0 = pos-npx_grd*(j_grid>0);
            pos1 = pos+npx_grd*(j_grid<npy_grd-1);
            
            
            CC[posC+(npx+j)*4]   = N0*fh[pos0] + N1*fh[pos] + N2*fh[pos1];
            CC[posC+(npx+j)*4+1] = N0*fhu[pos0] + N1*fhu[pos] + N2*fhu[pos1];
            CC[posC+(npx+j)*4+2] = N0*fhv[pos0] + N1*fhv[pos] + N2*fhv[pos1];
            CC[posC+(npx+j)*4+3] = H[j*npx+npx-1];
        }
        
    }


    delete [] fh;
    delete [] fhu;
    delete [] fhv;
}




// Devuelve 0 si todo ha ido bien, 1 si no existe algun fichero, y 2 si no hay memoria CPU suficiente
int cargarDatosProblema(string fich_ent, string &prefijo, double **dVh, double2 **dVhu, double **dVH, double **CC, double **tCC, unsigned long int &time_grid, double &xmin, double& xmax, double &ymin, double& ymax, int &npx, int &npy, int &nvol, double& dx, double& dy, double &tiempo_tot, double &tiempo_guardar, double &CFL, double& Q, double& DS, double& FR, double& Opar, double &OC, double &HQB, double &vmax, double &epsilon_h)
{
	int i, j;
    int pos;
    string file_bat;
    string file_eta;
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
	fich.close();
    
    //     SPL=5;  // Sponge layer
    FR/=8.0; //Chezy law
    HQB=5.0*DS;
    Q=DS*sqrt(G*(R-1.0)*DS)/(1.0-porosity);
    Opar=FR/((R-1.0)*G*DS);
    
#if Read_bathy==1
    
    // Reservamos memoria
    double *read_aux1;
    double *read_aux2;
    
    abrirbatiGRD((char *) file_bat.c_str());

    int npx_grd, npy_grd;
//     unsigned long int time_grid;

    leerTamanoMallaGRD((char *) file_bat.c_str(),&npx_grd, &npy_grd);
    double *xp;
    double *yp;
    xp=new double[npx_grd];
    yp=new double[npy_grd];
    leerLongitudGRD(xp);
    leerLatitudGRD(yp);
    double xmin_bathy, xmax_bathy, ymin_bathy, ymax_bathy;
    xmin_bathy=xp[0];
    xmax_bathy=xp[npx_grd-1];
    ymin_bathy=yp[0];
    ymax_bathy=yp[npy_grd-1];
    // Convertimos lat y longitud en distancia para sacar el dx,dy.
    //dy
    dy = dy_transform(yp[0],yp[1]);
    dx = dx_transform(xp[0],xp[1],yp[0]);
//     dy = 2*Re*(yp[1]-yp[0])/2.0 *PI/180;
//     dx = 2*Re*asin( sqrt( cos(yp[0]*PI/180)*cos(yp[0]*PI/180)*sin(0.5*(xp[1]-xp[0])*PI/180) *sin(0.5*(xp[1]-xp[0])*PI/180) )  );

    delete[] xp;
    delete[] yp;

    npx+=2*SPL;
    npy+=2*SPL;
    xmin_bathy-=SPL*(xmax_bathy-xmin_bathy)/(npx);
    ymin_bathy-=SPL*(ymax_bathy-ymin_bathy)/(npy);
    xmax_bathy+=SPL*(xmax_bathy-xmin_bathy)/(npx);
    ymax_bathy+=SPL*(ymax_bathy-ymin_bathy)/(npy);;
    
#if interpolate==0
    npx=npx_grd;
    npy=npy_grd;
#endif
    
    if ((npx != npx_grd || npy!= npy_grd) && (interpolate==0 && (Read_bathy==1 || Read_eta==1))){
        printf("You have set not to interpolate, but the number of points in data does not match the number of points in the reading file. Either set interpolation to 1 or set a proper number of points in data.dat. End of program");
        return 2;
    }

    
    // Descomentar esto
    xmin = 0.0;
    xmax = xmin + npx*dx;
    ymin = 0.0;
    ymax = ymin + npy*dy;

    nvol = npx*npy;


    
  //Quitar esto  
//     dx=(xmax-xmin)/(npx);
//     dy=(ymax-ymin)/(npy);
//     
//     xmin_bathy=xmin;
//     xmax_bathy=xmax;
//     ymin_bathy=ymin;
//     ymax_bathy=ymax;

//     
    int nvol_grd=npx_grd*npy_grd;

    read_aux1=new double[nvol_grd];
    read_aux2=new double[nvol_grd];
    
    double *inter_aux1, *inter_aux2;
    inter_aux1=new double[nvol]();
    inter_aux2=new double[nvol]();
    

    *dVhu = new double2 [nvol]();
    *dVh = new double [nvol]();
    *dVH = new double [nvol]();
        
	if (*dVH == NULL) {
		cerr << "Error: No hay memoria CPU suficiente" << endl;
        if (read_aux1 != NULL) delete[] read_aux1;
        if (read_aux2 != NULL) delete[] read_aux2;
        if (inter_aux1 != NULL) delete[] inter_aux1;
        if (inter_aux2 != NULL) delete[] inter_aux2;
 		liberarMemoria(*dVh,*dVhu,*dVH);
		return 2;
	}
	

    leerBatimetriaGRD(npx_grd, npy_grd, read_aux1);
#if Read_eta==1
    abrirCI_GRD(file_eta.c_str());
    leerCI_GRD(read_aux2, time_grid);
    int npx_grd_CI, npy_grd_CI;
    leerTamanoMallaGRD((char *) file_eta.c_str(),&npx_grd_CI, &npy_grd_CI);
    double *xp_CI;
    double *yp_CI;
    xp_CI=new double[npx_grd_CI];
    yp_CI=new double[npy_grd_CI];
    leerLongitudGRD(xp_CI,1);
    leerLatitudGRD(yp_CI,1);
    double xmin_CI, xmax_CI, ymin_CI, ymax_CI;
    xmin_CI=xp_CI[0];
    xmax_CI=xp_CI[npx_grd_CI-1];
    ymin_CI=yp_CI[0];
    ymax_CI=yp_CI[npy_grd_CI-1];

    
    delete[] xp_CI;
    delete[] yp_CI;
    
    
    
    *tCC=new double [time_grid]();
    *CC = new double [(2*npx+2*npy)*4*time_grid]();
    if (*CC == NULL) {
        cerr << "Error: No hay memoria CPU suficiente" << endl;
        liberarMemoria(*dVh,*dVhu,*dVH);
        return 2;
    }
    
#endif    
#if interpolate==1
    interpolate_fun(read_aux1,inter_aux1,npx,npy,npx_grd,npy_grd,xmin_bathy,xmax_bathy,ymin_bathy,ymax_bathy);
#endif
#if (Read_eta==1 && interpolate==1)
    interpolate_CI(read_aux2,inter_aux2,npx,npy,npx_grd_CI,npy_grd_CI,xmin_bathy,xmax_bathy,ymin_bathy,ymax_bathy,xmin_CI,xmax_CI,ymin_CI,ymax_CI, npx_grd_CI, npy_grd_CI);
#endif
#if Read_eta==1
    interpolate_CC(file_eta, (*CC), (*tCC), npx, npy, xmin_bathy, xmax_bathy, ymin_bathy, ymax_bathy, time_grid, xmin_CI,xmax_CI,ymin_CI,ymax_CI, inter_aux1);
    cerrarGRD(ncid_CI);
#endif
cerrarGRD(ncid_bati);

    for (j=0; j<npy; j++) {
        for (i=0; i<npx; i++) {
            pos=j*npx+i;
#if interpolate==1
            (*dVH)[pos]=(double) (inter_aux1[pos]) ;
#else
            (*dVH)[pos]=(double) (read_aux1[pos]);
#endif
//              printf("%d  %d   %.15f \n",i,j,(*dVH)[pos]);
#if Read_eta==0   
            double x,y;
            x = xmin + 0.5*dx + i*dx;
            y = ymin + 0.5*dy + j*dy;
            (*dVh)[pos] = h_inicial(x,y,(*dVH)[pos]) + (*dVH)[pos];
            velocidad(x,y,(*dVhu)[pos]);
            (*dVhu)[pos].x *= (*dVh)[pos];
            (*dVhu)[pos].y *= (*dVh)[pos];
#else
#if interpolate==0
            (*dVh)[pos] = (double) (read_aux2[pos]) *UNITS + (*dVH)[pos];
#else
            (*dVh)[pos] = (double) (inter_aux2[pos]) *UNITS + (*dVH)[pos];
//             (*dVh)[pos] =1.0;
#endif
//                          printf("%d  %d   %.15f \n",i,j,inter_aux2[pos]);
            (*dVh)[pos] *= ((*dVh)[pos]>=0.0);
            (*dVhu)[pos].x = 0.0;
            (*dVhu)[pos].y = 0.0;            
#endif
        }         
    }
    delete[] read_aux1;
    delete[] read_aux2;
    delete[] inter_aux1;
    delete[] inter_aux2;
 
    
// return 2;
    
#else  // Del Read_bathy==1
    
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

    *dVhu = new double2 [nvol]();
    *dVh = new double [nvol]();
    *dVH = new double [nvol]();
    
	if (*dVH == NULL) {
		cerr << "Error: No hay memoria CPU suficiente" << endl;
 		liberarMemoria(*dVh,*dVhu,*dVH);
		return 2;
	}


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
        }         
    }


#endif

	return 0;
}


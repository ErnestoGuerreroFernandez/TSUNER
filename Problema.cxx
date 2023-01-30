#include <sys/stat.h> 
#include <fstream>
#include <netcdf.h>
#include "cini.cxx"

int ncid_bati;
int ncid_emisores;

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
void abriremisoresGRD(const char *nombre_fich)
{
	int iret;

	iret = nc_open(nombre_fich, NC_NOWRITE, &ncid_emisores);
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
void leerLongitudGRD(double *lon)
{
	int id, iret;
	int lon_id;
	id = ncid_bati;
	// El nombre de la variable puede ser lon o x
	iret = nc_inq_varid(id, "lon", &lon_id);
	if (iret != NC_NOERR) {
		iret = nc_inq_varid(id, "x", &lon_id);
		check_errGRD(iret);
	}
	iret = nc_get_var_double(id, lon_id, lon);
	check_errGRD(iret);
}

void leerLatitudGRD(double *lat)
{
	int id, iret;
	int lat_id;

	id = ncid_bati;
	// El nombre de la variable puede ser lat o y
	iret = nc_inq_varid(id, "lat", &lat_id);
	if (iret != NC_NOERR) {
		iret = nc_inq_varid(id, "y", &lat_id);
		check_errGRD(iret);
	}
	iret = nc_get_var_double(id, lat_id, lat);
	check_errGRD(iret);
}
void leerBatimetriaGRD(int num_volx, int num_voly, float *bati)
{
	int id, iret;
	int z_id;
	const size_t start[2] = {0, 0};
	const size_t count[2] = {num_voly, num_volx};

	id = ncid_bati;
	iret = nc_inq_varid(id, "z", &z_id);
	check_errGRD(iret);
	iret = nc_get_vara_float(id, z_id, start, count, bati);
	check_errGRD(iret);
}

void leerEmisoresGRD(int num_volx, int num_voly, float *emisores)
{
	int id, iret;
	int z_id;
	const size_t start[2] = {0, 0};
	const size_t count[2] = {num_voly, num_volx};

	id = ncid_emisores;
	iret = nc_inq_varid(id, "z", &z_id);
	check_errGRD(iret);
	iret = nc_get_vara_float(id, z_id, start, count, emisores);
	check_errGRD(iret);
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

/*
// Devuelve 0 si todo ha ido bien, 1 si no existe algun fichero, y 2 si no hay memoria CPU suficiente
int cargarDatosProblema(string fich_ent, string &prefijo, double **dVh, double2 **dVhu,double **dVH, double &xmin, double& xmax, double &ymin, double& ymax, int &npx, int &npy, int &nvol, double& dx, double& dy, double &tiempo_tot, double &tiempo_guardar, double& dtmax, double &CFL, double& tg1, double& tg2, double& tg3, double& ca, double &vmax, double &epsilon_h, double& H, double& L, double& U, double& T)
{
	int i, j;
    double x,y;
    int pos;
    double prof;
    string file_bat;
    string file_emisores;
 	
	// Directorio donde se encuentran los ficheros de datos
	string directorio;

	// Ponemos en directorio el directorio donde est�n los ficheros de datos
	i = fich_ent.find_last_of("/");
	if (i > -1) {
		// Directorios indicados con '/' (S.O. distinto de windows)
		directorio = fich_ent.substr(0,i)+"/";
	}
	else {
		i = fich_ent.find_last_of("\\");
		if (i > -1) {
			// Directorios indicados con '\' (S.O. Windows)
			directorio = fich_ent.substr(0,i)+"\\";
		}
		else {
			// No se ha especificado ning�n directorio para los ficheros de datos
			directorio = "";
		}
	}
    
	// Leemos los datos del problema del fichero de entrada
	ifstream fich(fich_ent.c_str());

    
    fich>> file_bat;
    fich>> file_emisores;
	fich >> tiempo_tot;
	fich >> tiempo_guardar;
    fich >> dtmax;
	fich >> CFL;
	fich >> epsilon_h;
    fich >> tg1;
    fich >> tg2;
    fich >> tg3;
    fich >> ca;
	fich >> vmax;
    fich >> H;
    fich >> L;
    fich >> prefijo;
    U=sqrt(H*9.81);
    T=L/U;
	fich.close();
    
    tg1=tan(tg1*PI/180.0);
    tg1*=(L/H);
    tg2=tan(tg2*PI/180.0);
    tg2*=(L/H);
    tg3=tan(tg3*PI/180.0);
    tg3*=(L/H);

    epsilon_h/=H;

    vmax/=U;
    
    ca*=(L/H);
    // Reservamos memoria
    float *read_aux1;
    float *read_aux2;
    
    abrirbatiGRD((char *) file_bat.c_str());
    abriremisoresGRD((char *) file_emisores.c_str());
    leerTamanoMallaGRD((char *) file_bat.c_str(),&npx, &npy);
    double *xp;
    double *yp;
    xp=new double[npx];
    yp=new double[npy];
    leerLongitudGRD(xp);
    leerLatitudGRD(yp);
    xmin=xp[0];
    xmax=xp[npx-1];
    ymin=yp[0];
    ymax=yp[npy-1];
    delete[] xp;
    delete[] yp;


    nvol = npx*npy;

    dx=(xmax-xmin)/(npx-1);
    dy=(ymax-ymin)/(npy-1);
    xmin=xmin-0.5*dx;
    xmax=xmax+0.5*dx;
    ymin=ymin-0.5*dy;
    ymax=ymax+0.5*dy;

    read_aux1=new float[nvol];
    read_aux2=new float[nvol];

    *dVhu = (double2 *) malloc(nvol*sizeof(double2));
    *dVh = (double *) malloc(nvol*sizeof(double));
    *dVH = (double *) malloc(nvol*sizeof(double));
    *emisores=(int *) malloc(nvol*sizeof(int));
    
    
	if (*dVH == NULL) {
		cerr << "Error: No hay memoria CPU suficiente" << endl;
        if (read_aux1 != NULL) delete[] read_aux1;
        if (read_aux2 != NULL) delete[] read_aux2;
 		liberarMemoria(*emisores,*dVh,*dVhu,*dVH);
		return 2;
	}

    leerBatimetriaGRD(npx, npy, read_aux1);
    leerEmisoresGRD(npx, npy, read_aux2);
    
    cerrarGRD(ncid_emisores);
    cerrarGRD(ncid_bati);
    for (j=0; j<npy; j++) {
        for (i=0; i<npx; i++) {
            pos=j*npx+i;
            (*dVh)[pos] = 0.0;
            (*emisores)[pos]=(int) read_aux2[pos];
            (*dVH)[pos]=(double) (-read_aux1[pos]);
            (*dVhu)[pos].x = 0.0;
            (*dVhu)[pos].y = 0.0;
        }         
    }
    delete[] read_aux1;
    delete[] read_aux2;

    xmin/=L;
    xmax/=L;
    ymin/=L;
    ymax/=L;
    dx/=L;
    dy/=L;
    tiempo_tot/=T;
    tiempo_guardar/=T;
    dtmax/=T;
	return 0;
}

*/

// Devuelve 0 si todo ha ido bien, 1 si no existe algun fichero, y 2 si no hay memoria CPU suficiente
int cargarDatosProblema(string fich_ent, string &prefijo, double **dVh, double2 **dVhu,double **dVH, double &xmin, double& xmax, double &ymin, double& ymax, int &npx, int &npy, int &nvol, double& dx, double& dy, double &tiempo_tot, double &tiempo_guardar, double &CFL, double& Q, double& DS, double& FR, double& Opar, double &OC, double &HQB, int SPL, double &vmax, double &epsilon_h)
{
	int i, j;
    double x,y;
    int pos;
    double prof;
    string file_bat;
    string file_emisores;
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

	fich.close();
    
    dx=(xmax-xmin)/(npx);
    dy=(ymax-ymin)/(npy);
    
//     SPL=5;  // Sponge layer
    FR/=8.0; //Chezy law
    HQB=5.0*DS;
    Q=DS*sqrt(G*(R-1.0)*DS)/(1.0-porosity);
    Opar=FR/((R-1.0)*G*DS);
    
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

	return 0;
}


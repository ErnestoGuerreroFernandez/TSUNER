#include <string.h>
#include "Constantes.hxx"
#include "Problema.cxx"

/*****************/
/* Funciones GPU */
/*****************/

extern "C" int comprobarSoporteCUDA(int cudaDevice);
extern "C" void inicializarGPU();

extern "C" int Exner(double *dVh, double2 *dVhu, double *dVH, int *fixH, double xmin, double xmax, double ymin, double ymax, int npx, int npy, int nvol, double dx, double dy, char *prefijo,  double tiempo_tot, double tiempo_guardar, double CFL, double Q, double DS, double FR, double Opar, double OC, double HQB, double vmax, double epsilon_h, double* d, double *CC, double *tCC, long unsigned int timeCC, double *dyrup, float *time_dyrup, double &tiempo, int CudaDevice, double *tides_bd, float *tides_time, int tidesNPT);

/*********************/
/* Fin funciones GPU */
/*********************/

int main(int argc, char *argv[])
{
	double2 *dVhu;  // hux huy
    double *dVh;  // h
    double* dVH; //batimetria
    double *CC;
    double *tCC;  // Boundary conditions.
    int *fixH; // Fixed bathymetry
    double *dyrup; // Dynamic rupture
    float *time_dyrup; // time of dynamic rupture
    long unsigned int timeCC;
    double *tides_bd;
    float *tides_time;
    int tidesNPT;
	int soporteCUDA, err;
    int cudaDevice;
	int npx, npy, nvol;
	double tiempo_tot, tiempo_guardar;
    double dx,dy;
	double CFL;
	double vmax;  // Los tg son parametros para la friccion de Coulomb. Se pueden sustituir por parametros para la friccion de Manning.
    double xmin, xmax, ymin, ymax;
	double epsilon_h;
	string fich_ent;
	string prefijo;      // Prefijo de los ficheros NetCDF
	double tiempo_gpu;
    double OC, Opar, Q, HQB, DS, FR;
//     int SPL = 5;

    // Fichero de datos
	fich_ent = argv[1];
    if (argc==3) {
        cudaDevice=atoi(argv[2]);
    }
    else {
        cudaDevice=0;
    }
	if (! existeFichero(fich_ent)) {
		cerr << "Error: Folder not found. '" << fich_ent << "'" << endl;
		return EXIT_FAILURE;
	}

	// Comprobamos si la tarjeta grafica soporta CUDA
	soporteCUDA = comprobarSoporteCUDA(cudaDevice);
    switch(soporteCUDA) {
        case 1:
            cerr << "Error: There is no GPU" << endl;
            return EXIT_FAILURE;
            break;
        case 2:
            cerr << "Error: There is no GPU that supports CUDA" << endl;
            return EXIT_FAILURE;
            break;
        case 3:
            cerr << "N. Device erroneo" << endl;
            return EXIT_FAILURE;
            break;
    }
	inicializarGPU();

    err = cargarDatosProblema(fich_ent, prefijo, &dVh, &dVhu, &dVH, &CC, &tCC, &fixH, &dyrup, &time_dyrup, timeCC, xmin, xmax, ymin, ymax, npx, npy,nvol, dx, dy, tiempo_tot, tiempo_guardar, CFL, Q, DS, FR, Opar, OC, HQB, vmax, epsilon_h, &tides_bd, &tides_time, tidesNPT);
	if (err > 0)
		return EXIT_FAILURE;
	cout << scientific;
    
    // Sponge layer dumping function and values at the boundary
    double *d;
    
#if SPL>0
    d = new double [SPL]();
        for (int i=0; i<SPL; i++){
            d[i] = (double(SPL-i)/double(SPL))*(double(SPL-i)/double(SPL));
//             d[i] = 1.0;
        }
#endif


	err = Exner(dVh, dVhu, dVH, fixH, xmin, xmax, ymin, ymax, npx, npy, nvol, dx, dy, (char *) prefijo.c_str(), tiempo_tot, tiempo_guardar, CFL, Q,DS,FR,Opar, OC, HQB, vmax,epsilon_h, d, CC, tCC, timeCC, dyrup, time_dyrup, tiempo_gpu,cudaDevice, tides_bd, tides_time, tidesNPT);
	if (err > 0) {
		if (err == 1)
			cerr << "Error: Not enough GPU memory" << endl;
		else if (err == 2)
			cerr << "Error: Not enough CPU memory" << endl;
		liberarMemoria(dVh,dVhu,dVH,fixH, CC, tCC, dyrup, time_dyrup);
		return EXIT_FAILURE;
	}
	cout << endl << "Total simulation time: " << tiempo_gpu << " sec." << endl;
	liberarMemoria(dVh, dVhu, dVH,fixH,CC, tCC, dyrup, time_dyrup);
#if SPL>0
    delete [] d;
#endif
    cout<<"Normal End."<<endl;

	return EXIT_SUCCESS;
}


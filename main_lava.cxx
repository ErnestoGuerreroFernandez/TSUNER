#include <string.h>
#include "Constantes.hxx"
#include "Problema.cxx"

/*****************/
/* Funciones GPU */
/*****************/

extern "C" int comprobarSoporteCUDA(int cudaDevice);
extern "C" void inicializarGPU();

extern "C" int Lava(double *dVh, double2 *dVhu, double *dVH, double xmin, double xmax, double ymin, double ymax, int npx, int npy, int nvol, double dx, double dy, char *prefijo,  double tiempo_tot, double tiempo_guardar, double CFL, double Q, double DS, double FR, double Opar, double OC, double HQB, double vmax, double epsilon_h, double* d, double hl, double Hl, double2 ql, double hr, double Hr, double2 qr, double ht, double Ht, double2 qt, double hb, double Hb, double2 qb, double &tiempo, int CudaDevice);

/*********************/
/* Fin funciones GPU */
/*********************/

int main(int argc, char *argv[])
{
	double2 *dVhu;  // hux huy
    double *dVh;  // h
    double* dVH; //batimetria
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
		cerr << "Error: No se ha encontrado el fichero '" << fich_ent << "'" << endl;
		return EXIT_FAILURE;
	}

	// Comprobamos si la tarjeta grafica soporta CUDA
	soporteCUDA = comprobarSoporteCUDA(cudaDevice);
    switch(soporteCUDA) {
        case 1:
            cerr << "Error: No hay tarjeta grafica" << endl;
            return EXIT_FAILURE;
            break;
        case 2:
            cerr << "Error: No hay ninguna tarjeta grafica que soporte CUDA" << endl;
            return EXIT_FAILURE;
            break;
        case 3:
            cerr << "N. Device erroneo" << endl;
            return EXIT_FAILURE;
            break;
    }
	inicializarGPU();
    
    err = cargarDatosProblema(fich_ent, prefijo, &dVh, &dVhu, &dVH, xmin, xmax, ymin, ymax, npx, npy,nvol, dx, dy, tiempo_tot, tiempo_guardar, CFL, Q, DS, FR, Opar, OC, HQB, vmax, epsilon_h);
	if (err > 0)
		return EXIT_FAILURE;
	cout << scientific;
    
            // Sponge layer dumping function and values at the boundary
        double *d;
        double hl, Hl, hr, Hr, ht, Ht, hb, Hb;
        double2 ql, qr, qb, qt;
        
#if SPL>0
        d = new double [SPL]();
        
        int pmx, pmy;
        pmx = int(npx/2);
        pmy = int(npy/2);

        
        hl = dVh[pmy*npx];
        Hl = dVH[pmy*npx];
        ql.x = 0.0;
        ql.y = 0.0;

        hr = dVh[pmy*npx + npx-1];
        Hr = dVH[pmy*npx + npx-1];
        qr.x = 0.0;
        qr.y = 0.0;

        hb = dVh[pmx];
        Hb = dVH[pmx];
        qb.x = 0.0;
        qb.y = 0.0;

        ht = dVh[npx*npy-pmx];
        Ht = dVH[npx*npy-pmx];
        qt.x = 0.0;
        qt.y = 0.0;
        
        hl = 1.0;
        Hl = 0.0;
        hr = hl;
        Hr = Hl;
//         hb = 00;
//         Hb = Hl;
        ht = hl;
        Ht = Hl;

        for (int i=0; i<SPL; i++){
            d[i] = (double(SPL-i)/double(SPL))*(double(SPL-i)/double(SPL));
//             d[i] = 0.0;
        }

        
#endif

        
	err = Lava(dVh, dVhu, dVH, xmin, xmax, ymin, ymax, npx, npy, nvol, dx, dy, (char *) prefijo.c_str(), tiempo_tot, tiempo_guardar, CFL, Q,DS,FR,Opar, OC, HQB, vmax,epsilon_h, d, hl, Hl, ql, hr, Hr, qr, ht, Ht, qt, hb, Hb, qb, tiempo_gpu,cudaDevice);
	if (err > 0) {
		if (err == 1)
			cerr << "Error: No hay memoria GPU suficiente" << endl;
		else if (err == 2)
			cerr << "Error: No hay memoria CPU suficiente" << endl;
		liberarMemoria(dVh,dVhu,dVH);
		return EXIT_FAILURE;
	}
	cout << endl << "Tiempo: " << tiempo_gpu << " seg." << endl;
	liberarMemoria(dVh, dVhu, dVH);
    
    cout<<"Fin Normal."<<endl;

	return EXIT_SUCCESS;
}


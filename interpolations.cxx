#include <netcdf.h>
#include "Constantes.hxx" 


bool ErrorNetCDF_inte;

int AUX = 1;


void check_errGRD_inte(int iret)
{
        if ((iret != NC_NOERR) && (! ErrorNetCDF_inte)) {
                fprintf(stderr, "%s\n", nc_strerror(iret));
                ErrorNetCDF_inte = true;
        }
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

inline void dy_transform(double &val){
    val = (2*PI*Re)/360;
}

inline void dx_transform(double ymin, double dy, int tam, double xmax, double xmin, double *dx, double dx_cartesianas, double &dmax){
    dmax = -1.0;
    for (int i=0;i<tam;i++){
        double lat;
        lat = ymin + 0.5*dy + i*dy; 
        dx[i] = (2*PI*Re*cos(lat*PI/180))/360;
        dmax = max(dmax,dx[i]*(xmax-xmin));
        dx[i] *= dx_cartesianas;
    }
//     return 2*Re*asin( sqrt( cos(phi*PI/180)*cos(phi*PI/180)*sin(0.5*(alpha2-alpha1)*PI/180) *sin(0.5*(alpha2-alpha1)*PI/180) )  );
}


/*
void interpolate_bathy(double *bathy, double *p, int npx, int npy, int npx_grd, int npy_grd, double xmin, double xmax, double ymin, double ymax, double *dx_grid, double dy_grid){
    
    double x,y;
    double dx,dy;
//     double dx_grid, dy_grid;
    double xt, yt;
    int i_grid, j_grid;
    double N1,N2,N3,N4,N5,N6,N7,N8,N9;
    
    dx = (xmax-xmin)/npx;
    dy = (ymax-ymin)/npy;
    

//     dx_grid = (xmax-xmin)/npx_grd;
//     dy_grid = (ymax-ymin)/npy_grd;
    
    double dmax = xmax;
    

    for (int j=0; j<npy; j++){
        for (int i=0; i<npx; i++){
            int pos, posp;
            int posL, posR, posU, posD, posLU, posLD, posRU, posRD;
            
            x = xmin + 0.5*dx + i*dx;
            y = ymin + 0.5*dy + j*dy;
        
            xt = 1.0/(xmax-xmin)*(2*x-(xmax+xmin));
            yt = 1.0/(ymax-ymin)*(2*y-(ymax+ymin));
        
            double xmin_grid;

            j_grid = int((y - ymin - 0.5*dy_grid)/dy_grid);  // El limite de esto cuando dy tiende a cero es -0.5, que se redondea a 0. Por tanto, esto no dara enteros negativos.
            if(j_grid>npy_grd-1) j_grid=npy_grd-1;
            
            xmin_grid = (dmax - dx_grid[j_grid]*(npx_grd-1))/2;
            
            i_grid = int((x - xmin_grid - 0.5*dx_grid[j_grid])/dx_grid[j_grid]);
            if(i_grid>npx_grd-1) i_grid=npx_grd-1;
            
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


void interpolate_CI(double *u, double *p, int npx, int npy, int npx_grd, int npy_grd, double xmin, double xmax, double ymin, double ymax, int npx_grd_CI, int npy_grd_CI, double *dx_CI, double dy_CI, double *H){
    
    double x,y;
    double dx,dy;
    double dx_grid, dy_grid;
    double xt, yt;
    int i_grid, j_grid;
    double N1,N2,N3,N4,N5,N6,N7,N8,N9;
    
    dx = (xmax-xmin)/npx;
    dy = (ymax-ymin)/npy;
    
//     dx_grid = (xmax_CI-xmin_CI)/npx_grd_CI;
//     dy_grid = (ymax_CI-ymin_CI)/npy_grd_CI;
    
    double dmax = xmax;
    
    for (int j=0; j<npy; j++){
        for (int i=0; i<npx; i++){
            int pos, posp;
            int posL, posR, posU, posD, posLU, posLD, posRU, posRD;
            
            x = xmin + 0.5*dx + i*dx;
            y = ymin + 0.5*dy + j*dy;
        
            xt = 1.0/(xmax-xmin)*(2*x-(xmax+xmin));
            yt = 1.0/(ymax-ymin)*(2*y-(ymax+ymin));
        
            double xmin_CI;

            j_grid = int((y - ymin - 0.5*dy_CI)/dy_CI);
            if(j_grid>npy_grd-1) j_grid=npy_grd-1;
            
            xmin_CI = (dmax - dx_CI[j_grid]*(npx_grd-1))/2;
            
            i_grid = int((x - xmin_CI - 0.5*dx_CI[j_grid])/dx_CI[j_grid]);
            if(i_grid>npx_grd-1) i_grid=npx_grd-1;
            
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


        
            p[posp] = u[posLD]*N1 + u[posD]*N2 + u[posRD]*N3 + u[posL]*N4 + u[pos]*N5 + u[posR]*N6 + u[posLU]*N7 + u[posU]*N8 + u[posRU]*N9; 
            
            if(p[posp]<-1e20) p[posp] = H[posp];
            if(p[posp]>1e20) p[posp] = H[posp];
            
        }    
    }
}
*/

void interpolate_CC(string file, double *CC,  double *tCC, int npx, int npy, double xmin, double xmax, double ymin, double ymax, double xmin_grid, double ymin_grid, long unsigned int time, double dx_grid, double dy_grid, double *H, int id){
    
    
    double x,y;
    double dx,dy;
//     double dx_grid, dy_grid;
    double xt, yt;
    int i_grid, j_grid;
    int idu;
    int var_id, dimId;
// 	int id = ncid_CI;
    long unsigned int npx_grd, npy_grd;
    

    
    double *fh, *fhu, *fhv;
    

	
    // We build the file XXXu.nc
    file.erase(file.end()-4);
    file.insert(file.end()-3,'u');
	
	int iret = nc_open(file.c_str(), NC_NOWRITE, &idu);

    iret = nc_inq_dimid(idu, "lat", &dimId);
        if (iret != NC_NOERR) {
        iret = nc_inq_dimid(idu, "LAT", &dimId);
		check_errGRD_inte(iret);
    }
    iret = nc_inq_dimlen(idu, dimId, &npy_grd);
    check_errGRD_inte(iret);
    
    iret = nc_inq_dimid(idu, "lon", &dimId);
    if (iret != NC_NOERR) {
        iret = nc_inq_dimid(idu, "LON", &dimId);
		check_errGRD_inte(iret);
    }
    iret = nc_inq_dimlen(idu, dimId, &npx_grd);
    check_errGRD_inte(iret);
    
    iret = nc_inq_dimid(idu, "time", &dimId);
    if (iret != NC_NOERR) {
        iret = nc_inq_dimid(idu, "TIME", &dimId);
		check_errGRD_inte(iret);
    }
    check_errGRD_inte(iret);
    
    static size_t start_t[1] = {0};
	static size_t count_t[1] = {time};
    
    iret = nc_inq_varid(idu, "time", &var_id);
    if (iret != NC_NOERR) {
        iret = nc_inq_varid(idu, "TIME", &var_id);
		check_errGRD_inte(iret);
    }
	iret = nc_get_vara_double(idu, var_id, start_t, count_t, tCC);
	check_errGRD_inte(iret);
    
    
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
    if (iret != NC_NOERR) {
        iret = nc_inq_varid(idu, "UA", &var_id);
		check_errGRD_inte(iret);
    }
	iret = nc_get_vara_double(idu, var_id, start, count, fhu);
	check_errGRD_inte(iret);
    iret = nc_close(idu);
	check_errGRD_inte(iret);
    
    // We build the file XXXv.nc
    file.erase(file.end()-4);
    file.insert(file.end()-3,'v');
	
	iret = nc_open(file.c_str(), NC_NOWRITE, &idu);
	check_errGRD_inte(iret);
    iret = nc_inq_varid(idu, "va_lc", &var_id);
    if (iret != NC_NOERR) {
        iret = nc_inq_varid(idu, "VA", &var_id);
		check_errGRD_inte(iret);
    }
	iret = nc_get_vara_double(idu, var_id, start, count, fhv);
	check_errGRD_inte(iret);
    iret = nc_close(idu);
	check_errGRD_inte(iret);

    
	iret = nc_inq_varid(id, "ha_lc", &var_id);
    if (iret != NC_NOERR) {
	iret = nc_inq_varid(id, "HA", &var_id);
		check_errGRD_inte(iret);
    }
	iret = nc_get_vara_double(id, var_id, start, count, fh);
	check_errGRD_inte(iret);


    dx = (xmax-xmin)/npx;
    dy = (ymax-ymin)/npy;
    

    
    double dmax = xmax;
    
    for (int t=0; t<time; t++){
        int pos, pos0, pos1;
        int posC = t*(2*npx+2*npy)*4;
        double N0, N1, N2;
        for (int i=0; i<npx; i++){
            x = xmin + 0.5*dx + i*dx;
            y = ymin ;
            xt = 1.0/(xmax-xmin)*(2*x-(xmax+xmin));
            
            N0 = L0(xt);
            N1 = L1(xt);
            N2 = L2(xt);
            
            i_grid = int((x - xmin_grid - 0.5*dx_grid)/dx_grid);
            if(i_grid>npx_grd-1) i_grid=npx_grd-1;
            j_grid = int((y - ymin_grid - 0.5*dy_grid)/dy_grid);
            if(j_grid>npy_grd-1) j_grid=npy_grd-1;
            
            
            
            // Frontera inf
            pos = t*npx_grd*npy_grd + npx_grd*j_grid + i_grid;
            pos0 = pos-1*(i_grid>0);
            pos1 = pos+1*(i_grid<npx_grd-1);
            
            CC[posC+i*4] = N0 * fh[pos0] + N1*fh[pos] + N2*fh[pos1];
            CC[posC+i*4+1] = N0 * fhu[pos0] + N1*fhu[pos] + N2*fhu[pos1];
            CC[posC+i*4+2] = N0 * fhv[pos0] + N1*fhv[pos] + N2*fhv[pos1];
            CC[posC+i*4+3] = H[i];
            
            if(CC[posC+i*4]<-1e20) CC[posC+i*4] = H[i]; 
            if(CC[posC+i*4+1]<-1e20) CC[posC+i*4+1] = 0.0; 
            if(CC[posC+i*4+2]<-1e20) CC[posC+i*4+2] = 0.0;
            if(CC[posC+i*4]>1e20) CC[posC+i*4] = H[i]; 
            if(CC[posC+i*4+1]>1e20) CC[posC+i*4+1] = 0.0; 
            if(CC[posC+i*4+2]>1e20) CC[posC+i*4+2] = 0.0;
            
            // Frontera sup
            y = ymax ;
            i_grid = int((x - xmin_grid - 0.5*dx_grid)/dx_grid);
            if(i_grid>npx_grd-1) i_grid=npx_grd-1;
            j_grid = int((y - ymin_grid - 0.5*dy_grid)/dy_grid);
            if(j_grid>npy_grd-1) j_grid=npy_grd-1;
            pos = t*npx_grd*npy_grd + (j_grid)*npx_grd + i_grid;
            pos0 = pos-1*(i_grid>0);
            pos1 = pos+1*(i_grid<npx_grd-1);
            
            
            CC[posC+(npx+npy+i)*4]   = N0*fh[pos0] + N1*fh[pos] + N2*fh[pos1];
            CC[posC+(npx+npy+i)*4+1] = N0*fhu[pos0] + N1*fhu[pos] + N2*fhu[pos1];
            CC[posC+(npx+npy+i)*4+2] = N0*fhv[pos0] + N1*fhv[pos] + N2*fhv[pos1];
            CC[posC+(npx+npy+i)*4+3] = H[(npy-1)*npx+i];
            
            if(CC[posC+(npx+npy+i)*4]<-1e20) CC[posC+(npx+npy+i)*4] = H[(npy-1)*npx+i]; 
            if(CC[posC+(npx+npy+i)*4+1]<-1e20) CC[posC+(npx+npy+i)*4+1] = 0.0; 
            if(CC[posC+(npx+npy+i)*4+2]<-1e20) CC[posC+(npx+npy+i)*4+2] = 0.0;
            if(CC[posC+(npx+npy+i)*4]>1e20) CC[posC+(npx+npy+i)*4] = H[(npy-1)*npx+i]; 
            if(CC[posC+(npx+npy+i)*4+1]>1e20) CC[posC+(npx+npy+i)*4+1] = 0.0; 
            if(CC[posC+(npx+npy+i)*4+2]>1e20) CC[posC+(npx+npy+i)*4+2] = 0.0;
            
        }
        for (int j=0; j<npy; j++){
            y = ymin + 0.5*dy + j*dy;
            x = xmin ;
            yt = 1.0/(ymax-ymin)*(2*y-(ymax+ymin));
            
            N0 = L0(yt);
            N1 = L1(yt);
            N2 = L2(yt);
            
            i_grid = int((x - xmin_grid - 0.5*dx_grid)/dx_grid);
            if(i_grid>npx_grd-1) i_grid=npx_grd-1;
            j_grid = int((y - ymin_grid - 0.5*dy_grid)/dy_grid);
            if(j_grid>npy_grd-1) j_grid=npy_grd-1;
            
            
            // Frontera izq
            pos = t*npx_grd*npy_grd + j_grid*npx_grd + i_grid;
            pos0 = pos-npx_grd*(j_grid>0);
            pos1 = pos+npx_grd*(j_grid<npy_grd-1);
            
            CC[posC+(2*npx+npy+j)*4]   = N0*fh[pos0] + N1*fh[pos] + N2*fh[pos1];
            CC[posC+(2*npx+npy+j)*4+1] = N0*fhu[pos0] + N1*fhu[pos] + N2*fhu[pos1];
            CC[posC+(2*npx+npy+j)*4+2] = N0*fhv[pos0] + N1*fhv[pos] + N2*fhv[pos1];
            CC[posC+(2*npx+npy+j)*4+3] = H[j*npx];
            
            if(CC[posC+(2*npx+npy+j)*4]<-1e20) CC[posC+(2*npx+npy+j)*4] = H[j*npx]; 
            if(CC[posC+(2*npx+npy+j)*4+1]<-1e20) CC[posC+(2*npx+npy+j)*4+1] = 0.0; 
            if(CC[posC+(2*npx+npy+j)*4+2]<-1e20) CC[posC+(2*npx+npy+j)*4+2] = 0.0;
            if(CC[posC+(2*npx+npy+j)*4]>1e20) CC[posC+(2*npx+npy+j)*4] = H[j*npx]; 
            if(CC[posC+(2*npx+npy+j)*4+1]>1e20) CC[posC+(2*npx+npy+j)*4+1] = 0.0; 
            if(CC[posC+(2*npx+npy+j)*4+2]>1e20) CC[posC+(2*npx+npy+j)*4+2] = 0.0;
            
           
            
            // Frontera drch
            x = xmax ;
            i_grid = int((x - xmin_grid - 0.5*dx_grid)/dx_grid);
            if(i_grid>npx_grd-1) i_grid=npx_grd-1;
            pos = t*npx_grd*npy_grd + j_grid*npx_grd + i_grid;
            pos0 = pos-npx_grd*(j_grid>0);
            pos1 = pos+npx_grd*(j_grid<npy_grd-1);
            
            
            CC[posC+(npx+j)*4]   = N0*fh[pos0] + N1*fh[pos] + N2*fh[pos1];
            CC[posC+(npx+j)*4+1] = N0*fhu[pos0] + N1*fhu[pos] + N2*fhu[pos1];
            CC[posC+(npx+j)*4+2] = N0*fhv[pos0] + N1*fhv[pos] + N2*fhv[pos1];
            CC[posC+(npx+j)*4+3] = H[j*npx+npx-1];
            
            if(CC[posC+(npx+j)*4]<-1e20) CC[posC+(npx+j)*4] = H[j*npx+npx-1]; 
            if(CC[posC+(npx+j)*4+1]<-1e20) CC[posC+(npx+j)*4+1] = 0.0; 
            if(CC[posC+(npx+j)*4+2]<-1e20) CC[posC+(npx+j)*4+2] = 0.0;
            if(CC[posC+(npx+j)*4]>1e20) CC[posC+(npx+j)*4] = H[j*npx+npx-1]; 
            if(CC[posC+(npx+j)*4+1]>1e20) CC[posC+(npx+j)*4+1] = 0.0; 
            if(CC[posC+(npx+j)*4+2]>1e20) CC[posC+(npx+j)*4+2] = 0.0;
        }
        
    }
//     for (int t=0;t<time;t++){
//         int posC = t*(2*npx+2*npy)*4;
//         for (int j=0; j<npy; j++){
//              cout<<CC[posC+(2*npx+npy+j)*4]<<endl;
//         }
//     }
// exit(-1);

    delete [] fh;
    delete [] fhu;
    delete [] fhv;
}

void interpolate_spherics(double xmin, double xmax, double ymin, double ymax, double xmin_fun, double ymin_fun, int npx_fun, int npy_fun, double dx_fun, double dy_fun, int npx, int npy, double *fun, double *sol){
    // xmin, xmax, ... de sol
        // Interpolates fun into a mesh npx x npy,
    
    double x,y;
    double dx,dy;
//     double dx_fun, dy_fun;
    double xt, yt;
    int i_grid, j_grid;
    double N1,N2,N3,N4,N5,N6,N7,N8,N9;
    
    dx = (xmax-xmin)/(npx);
    dy = (ymax-ymin)/(npy);


    
    for (int j=0; j<npy; j++){
        for (int i=0; i<npx; i++){
            int pos, posp;
            int posL, posR, posU, posD, posLU, posLD, posRU, posRD;
            
            x = xmin + 0.5*dx + i*dx;
            y = ymin + 0.5*dy + j*dy;
        
            xt = 1.0/(xmax-xmin)*(2*x-(xmax+xmin));
            yt = 1.0/(ymax-ymin)*(2*y-(ymax+ymin));
        

            i_grid = int((x - xmin_fun - 0.5*dx_fun)/dx_fun);
            j_grid = int((y - ymin_fun - 0.5*dy_fun)/dy_fun);
            if(i_grid>npx_fun-1) i_grid=npx_fun-1;
            if(j_grid>npy_fun-1) j_grid=npy_fun-1;
            if(i_grid<0) i_grid=0;
            if(j_grid<0) j_grid=0;
            
            pos = j_grid*npx_fun + i_grid;
            posp = j*npx +i;
        
            posL = pos-1;
            posR = pos+1;
            posU = pos+npx_fun;
            posD = pos-npx_fun;
            posLU = posU-1;
            posLD = posD-1;
            posRU = posU+1;
            posRD = posD+1;
        
            if(i_grid==0){
                posL=pos;
                posLU=pos;
                posLD=pos;
            }
            if(i_grid==npx_fun-1){
                posR=pos;
                posRD=pos;
                posRU=pos;
            }
                
            if(j_grid==0){
                posD=pos;
                posLD=pos;
                posRD=pos;
            }
                
            if(j_grid==npy_fun-1){
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

        
            sol[posp] = fun[posLD]*N1 + fun[posD]*N2 + fun[posRD]*N3 + fun[posL]*N4 + fun[pos]*N5 + fun[posR]*N6 + fun[posLU]*N7 + fun[posU]*N8 + fun[posRU]*N9; 
            
        }    
    }
//     for (int j=0; j<npy; j++){
//         int pos = j*npx;
//         cout<<sol[pos]<<endl;
//     }
//     exit(-1);
}

void fix_eta(double *p, double *H, int npx, int npy){
    for (int j=0; j<npy; j++){
        for (int i=0; i<npx; i++){
            int posp = j*npx +i;
            if(p[posp]<-1e20) p[posp] = H[posp];
            if(p[posp]>1e20) p[posp] = H[posp];
        }
    }
}


void spheric2cartesian(double xmin, double xmax, double ymin, double ymax, double *dx, double dy, int npx, int npy, double *fun, double *sol){
    
    double dxG, dyG;
    double x,y;
    double xt, yt;
    int i_grid, j_grid;
    double N1,N2,N3,N4,N5,N6,N7,N8,N9;
    
    dxG = (xmax-xmin)/npx;
    dyG = (ymax-ymin)/npy;
    

    
    for (int i=0; i<npx; i++){
        for (int j=0; j<npy; j++){
            int pos, posp;
            int posL, posR, posU, posD, posLU, posLD, posRU, posRD;
            
            x = xmin + 0.5*dxG + i*dxG;
            y = ymin + 0.5*dyG + j*dyG;
        
            xt = 1.0/(xmax-xmin)*(2*x-(xmax+xmin));
            yt = 1.0/(ymax-ymin)*(2*y-(ymax+ymin));
        
            j_grid = int((y - ymin - 0.5*dy)/dy);
            if(j_grid>npy-1) j_grid=npy-1;
            if(j_grid<0) j_grid=0;
            
            double xminE = 0 + (xmax - dx[j_grid]*npx)/2;
            i_grid = int((x - xminE - 0.5*dx[j_grid])/dx[j_grid]);
            if(i_grid>npx-1) i_grid=npx-1;
            if(i_grid<0) i_grid=0;
            
            pos = j_grid*npx + i_grid;
            posp = j*npx +i;
            
            posL = pos-1;
            posR = pos+1;
            posU = pos+npx;
            posD = pos-npx;
            posLU = posU-1;
            posLD = posD-1;
            posRU = posU+1;
            posRD = posD+1;
        
            if(i_grid==0){
                posL=pos;
                posLU=pos;
                posLD=pos;
            }
            if(i_grid==npx-1){
                posR=pos;
                posRD=pos;
                posRU=pos;
            }
                
            if(j_grid==0){
                posD=pos;
                posLD=pos;
                posRD=pos;
            }
                
            if(j_grid==npy-1){
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
                
            
            sol[posp] = fun[posLD]*N1 + fun[posD]*N2 + fun[posRD]*N3 + fun[posL]*N4 + fun[pos]*N5 + fun[posR]*N6 + fun[posLU]*N7 + fun[posU]*N8 + fun[posRU]*N9; 
//             if(AUX==0) sol[posp] = fun[pos];
            
            
        }
    }
    
}


void spheric2cartesian_CC(double xmin, double xmax, double ymin, double ymax, double *dxE, double dyE, int npx, int npy, int time_grid, double *CC){
    
    
    double *aux;
    double x,y;
    double dx, dy;
    double N0, N1, N2;
    int pos, pos0, pos1, posA;
    double xt, yt;
    
    aux = new double[(2*npx+2*npy)*4]();
    
    dx = (xmax-xmin)/npx;
    dy = (ymax-ymin)/npy;

    for (int t=0; t<time_grid; t++){
        
        // Frontera sup e inf
        for(int i=0;i<npx;i++){
            int i_grid;
            // Frontera inf
            x = xmin + 0.5*dx + i*dx;
            y = ymin;
            
            double xminE = 0 + (xmax-dxE[0]*npx)/2;
            i_grid = int((x-xminE-0.5*dxE[0])/dxE[0]);
            if(i_grid>npx-1) i_grid=npx-1;
            if(i_grid<0) i_grid=0;
                        
            pos = t*(2*npx+2*npy)*4 + i_grid*4;
            pos0 = pos-4*(i_grid>0);
            pos1 = pos+4*(i_grid<npx-1);
                
            posA = i*4;
                        
            N0 = L0(xt);
            N1 = L1(xt);
            N2 = L2(xt);
            
            aux[posA]   = N0*CC[pos0] + N1*CC[pos] + N2*CC[pos1];
            aux[posA+1] = N0*CC[pos0+1] + N1*CC[pos+1] + N2*CC[pos1+1];
            aux[posA+2] = N0*CC[pos0+2] + N1*CC[pos+2] + N2*CC[pos1+2];
            aux[posA+3] = N0*CC[pos0+3] + N1*CC[pos+3] + N2*CC[pos1+3];
            
            
            // Frontera sup
            x = xmin + 0.5*dx + i*dx;
            y = ymax;
            
            xminE = 0 + (xmax-dxE[npy-1]*npx)/2;
            i_grid = int((x-xminE-0.5*dxE[npy-1])/dxE[npy-1]);
            if(i_grid>npx-1) i_grid=npx-1;
            if(i_grid<0) i_grid=0;
            
            
            pos = t*(2*npx+2*npy)*4 + (npx+npy)*4 + i_grid*4;
            pos0 = pos-4*(i_grid>0);
            pos1 = pos+4*(i_grid<npx-1);
                
            posA = (npx+npy)*4 + i*4;
                        
            N0 = L0(xt);
            N1 = L1(xt);
            N2 = L2(xt);
            
            aux[posA]   = N0*CC[pos0] + N1*CC[pos] + N2*CC[pos1];
            aux[posA+1] = N0*CC[pos0+1] + N1*CC[pos+1] + N2*CC[pos1+1];
            aux[posA+2] = N0*CC[pos0+2] + N1*CC[pos+2] + N2*CC[pos1+2];
            aux[posA+3] = N0*CC[pos0+3] + N1*CC[pos+3] + N2*CC[pos1+3];
        }
     
        // Frontera izq y drch
        for(int j=0;j<npy;j++){
            int j_grid;
            // Frontera izq
            x = xmin;
            y = ymin + 0.5*dy + j*dy;
            
            double yminE = 0 + (ymax-dyE*npy)/2;
            j_grid = int((y-yminE-0.5*dyE)/dyE);
            if(j_grid>npy-1) j_grid=npy-1;
            if(j_grid<0) j_grid=0;
                        
            pos = t*(2*npx+2*npy)*4 + (npx*2 + npy)*4 + j_grid*4;
            pos0 = pos-4*(j_grid>0);
            pos1 = pos+4*(j_grid<npy-1);
                
            posA = (npx*2+npy)*4 + j*4;
                        
            N0 = L0(xt);
            N1 = L1(xt);
            N2 = L2(xt);
            
            aux[posA]   = N0*CC[pos0] + N1*CC[pos] + N2*CC[pos1];
            aux[posA+1] = N0*CC[pos0+1] + N1*CC[pos+1] + N2*CC[pos1+1];
            aux[posA+2] = N0*CC[pos0+2] + N1*CC[pos+2] + N2*CC[pos1+2];
            aux[posA+3] = N0*CC[pos0+3] + N1*CC[pos+3] + N2*CC[pos1+3];
            
            
            // Frontera drch
            x = xmax;
            y = ymin + 0.5*dy + j*dy;
            
            yminE = 0 + (ymax-dyE*npy)/2;
            j_grid = int((y-yminE-0.5*dyE)/dyE);
            if(j_grid>npy-1) j_grid=npy-1;
            if(j_grid<0) j_grid=0;
            
            
            pos = t*(2*npx+2*npy)*4 + (npx)*4 + j_grid*4;
            pos0 = pos-4*(j_grid>0);
            pos1 = pos+4*(j_grid<npx-1);
                
            posA = (npx)*4 + j*4;
                        
            N0 = L0(xt);
            N1 = L1(xt);
            N2 = L2(xt);
            
            aux[posA]   = N0*CC[pos0] + N1*CC[pos] + N2*CC[pos1];
            aux[posA+1] = N0*CC[pos0+1] + N1*CC[pos+1] + N2*CC[pos1+1];
            aux[posA+2] = N0*CC[pos0+2] + N1*CC[pos+2] + N2*CC[pos1+2];
            aux[posA+3] = N0*CC[pos0+3] + N1*CC[pos+3] + N2*CC[pos1+3];
        }
        
//         if(AUX==1){
            for (int k=0;k<(npx*2+npy*2)*4;k++){
                CC[t*(npx*2+npy*2)*4+k] = aux[k];
            }
                
//         }
    }
    
    delete [] aux;
}














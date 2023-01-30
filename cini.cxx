#ifndef _cini_c_
#define _cini_c_

#include "Constantes.hxx"
#include <sys/stat.h>
#include <fstream>

double profundidad(double x, double y) {
    double H;
    double alfa=tan(15.0*PI/180);
    H = (1.4+x*alfa)*(x<0.0)+1.4*(x>=0.0);
    H = -10.0*exp(- 1.0* ((x+0)*(x+0) + (y+0)*(y+0)) );
    return H;
}


double h_inicial (double x, double y) {
    double h1;
    h1= profundidad(x,y)+1.0 + 15.0*exp(-5.0* ((x-2)*(x-2) + (y-2)*(y-2)));
    h1=h1*(h1>0.0);
    return h1;
}


void velocidad(double x , double y, double2& u) {
    u.x = x*x *0;
    u.y = y*y*y *0 ;
}


#endif


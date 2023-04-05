#ifndef _cini_c_
#define _cini_c_

#include "Constantes.hxx"
#include <sys/stat.h>
#include <fstream>

double profundidad(double x, double y) {
    double H;
    double alfa=tan(15.0*PI/180);
    H = (1.4+x*alfa)*(x<0.0)+1.4*(x>=0.0);
    H = -3.0*exp(- 5.0* ((x-0)*(x-0) + (y-0)*(y-0)) )  - -y*(y<-1.0) ;
//     H = -5.0*exp(- 5.0* ((x-2)*(x-2)*0 +  (y-2)*(y-2)) ) -5.0*exp(- 5.0* ((x+2)*(x+2) *0+ (y+2)*(y+2)) )  ; 
//     H = -3.0*exp(- 5.0* ((x-2)*(x-2) + (y-2)*(y-2)) ) -3.0*exp(- 5.0* ((x-2)*(x-2) + (y+2)*(y+2)) ) ;
//     H= -2.0*(x>1.0 && x<3.0);
//      H=-3.0*exp(- 5.0* ((x-0)*(x-0)  + (y-0)*(y-0)) )  - 0* 3.0*exp(- 5.0* ((x-0)*(x-0) + (y+2)*(y+2)) ) ;
//      H -= 3.0*exp(- 5.0* ((x-2)*(x-2) + (y-0)*(y-0)) ) + 3.0*exp(- 5.0* ((x+2)*(x+2) + (y-0)*(y-0)) );
//      H=0.0;
//      H = - x *(x>2.0);
    return H;
}


double h_inicial (double x, double y, double zb=0.0) {
    double h1;
    h1= profundidad(x,y) +1.0 +   1.0*exp(-5.0* ((x+2)*(x+2)   +   (y-3)*(y-3)))  + 0*35.0*exp(-5.0* ((x-0)*(x-0) + (y+2)*(y+2)));
//     h1=profundidad(x,y) + 10.0 + 755.0*exp(- 5.0* ((x-233)*(x-233)  +  (y-45)*(y-45)) )  ;
    
    
//     h1= profundidad(x,y) + 1.0 +  35.0*exp(-45.0* ((x-232)*(x-232)   + (y-44)*(y-44)))  + 0*35.0*exp(-5.0* ((x-0)*(x-0) + (y+2)*(y+2)));
    
    h1=h1*(h1>0.0);
    return h1;
}


void velocidad(double x , double y, double2& u) {
    u.x =   x*x*x;
    u.y = u.x *0 ;
    if(x<0) u.x = -u.x;
    u.x = x*0;
}


#endif


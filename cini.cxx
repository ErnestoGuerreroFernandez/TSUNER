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




#include "TSUNER.h"


double bathy(double x, double y) {  // Initial bathymetry function, H.
    double H;
    double alfa=tan(15.0*PI/180);

    H=-2.5*exp(- 5.0* ((x+0)*(x+0)  + (y-0)*(y-0)) )  - 0* 3.0*exp(- 5.0* ((x-0)*(x-0) + (y+2)*(y+2)) ) ;
    H += - x *(x>0.0);
    // H= 0.0;

    return H;
}


double initial_h (double x, double y) {  // Initial water depth, h.
    double h1;

    h1= bathy(x,y) +2.0 + 6*exp(-5.0* ((x+2)*(x+2)   +    (y-0)*(y-0)));
    // h1 = bathy(x,y) + 2.0;

    h1=h1*(h1>0.0);
    return h1;
}


void velocity(double x , double y, double2& u) {  // Initial velocity
    u.x = 0;
    u.y = 0;
}


int fixed_bed_CI(double x, double y){  // Set if the bathymetry moves (1) or not (0).
    int val;
    
    // if(x>=-4.0) val = 1;
    // else val=0;
    
    val=1;
    
    return val;
    
}




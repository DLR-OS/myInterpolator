//
// This file is part of myInterpolator - image void filling
// Copyright (C) 2022  Dirk 'jtk' Frommholz, DLR OS-SEC
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

#include "coolmath.h"


/****************************************************************************/


double common::rad2deg(double f_radAngle) {
    return f_radAngle*180/M_PI;
}


/****************************************************************************/


double common::deg2rad(double f_degAngle) {
    return f_degAngle*M_PI/180;
}


/****************************************************************************/


double common::gemanMcClure(double f_s, double f_x) {
    return square(f_x)/(square(f_s)+square(f_x));
}


/****************************************************************************/


double common::gaussianDistribution(double f_x, double f_mu, double f_sigma) {
    return (1/(f_sigma*sqrt(2*M_PI)))*exp(
        -(square(f_x-f_mu)/(2*square(f_sigma))
    ));                                  
}


/****************************************************************************/


double common::gaussianCDF(double f_x) {
    return 0.5*erfc(-f_x*0.7071067811865475244);
}


/****************************************************************************/

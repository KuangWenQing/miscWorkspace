#ifndef __OPENAPI_H__
#define __OPENAPI_H__

#include "myfunc.h"

/* get satellite position , velocity and DOPP
 * 
 * args   : time_t  time       I   second of week (s)
 *          eph_t   *eph       I   broadcast ephemeris
 *          double *usrpos     I   user position  {lat,lon,high} (rad | m)
 *          double *rs         O   sat position and velocity (ecef)
 *                                 {x,y,z,vx,vy,vz} (m|m/s)
 *          double *dopp       O   satellite DOPP (Hz)
 * 
 * return : none
 * *************************************************/
void satellite_pos_vel_dopp(time_t time, const eph_t *eph, const double *usrpos, double *rs, double *dopp);


/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite elevation angle
* args   : double *usrpos   I   user position, geodetic position {lat,lon,h} (rad,m)
*          double *satpos   I   satellite position, 
* 
* return : elevation angle (rad)  -pi/2<= el <=pi/2  is right
* ----------------------------------------------*/
double satellite_elevation(const double *usrpos, const double* satpos);


/* split the double type data to TYP_F32 type
 *  args    :   double      I   dData
 * 
 * 
 * 	return  :   A TYP_F32 type data
 * -------------------*/
TYP_F32 SplitExtF64(double dData);

#endif

#ifndef __OPENAPI_H__
#define __OPENAPI_H__

#include "myfunc.h"

/* get satellite position and DOPP
 * 
 * args   : time_t time        I   second of week (s)
 *          alm_t   *alm       I   almanac
 *          double *rs         O  sat position (ecef) and velocity {x,y,z, vx,vy,vz} (m | m/s)
 * 
 * return : none
 * *************************************************/
void satAlm_calc_satPos(time_t time, const alm_t *alm, double *rs);



/* compute satellite position by Ephemeris
 * 
 * args   : time_t time        I   second of week (s)
 *          eph_t   *eph       I   broadcast ephemeris
 *          double *rs         O  sat position (ecef) and velocity {x,y,z, vx,vy,vz} (m | m/s)
 * 
 * return : none
 * *************************************************/
void satEph_calc_satPos(time_t time, const eph_t *eph, double *rs);




/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite elevation angle and DOPP
* args   : double *usrpos   I   user position, geodetic position {lat,lon,h} (degree,m)
*          double *rs       I   satellite position, velocity {x,y,z, vx,vy,vz} (m, m/s)
*          double *dopp     IO  satellite DOPP (Hz)
* 
* return : elevation angle (rad)  -pi/2<= el <=pi/2  is right
* note :  user position unit  is  degree , elevation unit  is  rad
* ----------------------------------------------*/
double satellite_elevation_dopp(const double *usrpos, const double* rs, double *dopp);



/* split the double type data to TYP_F32 type
 *  args    :   double      I   dData
 * 
 * 	return  :   TYP_F32 type data
 * -------------------*/
TYP_F32 SplitExtF64(double dData);

#endif

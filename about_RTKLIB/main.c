#include <stdio.h>

#include "openAPI.h"

int main(int argc, char **argv)
{
	time_t time = 1112399999 % 604800 + 0.91728716015407952;

	eph_t test_eph = {
		3,
//		83,
//		595,
//		0,
//		0,
		1316,
//		1,
//		0,
		1112400000%604800,
//		1112400000%604800,
//      0,
		26560940.634528071,
		0.00673579110298,
		0.92743379988899999,
		0.53549319293800002,
		0.60389896875899995,
		2.4711168199300002,
		5.3766524565899999e-09,
		-8.27891621924e-09,
		-1.5250635476699999e-10,
		215.875,
		19.6875,
		1.0188668966299999e-06,
		7.5642019510299997e-06,
		-1.00582838058e-07,
		-6.51925802231e-08,
		518400,			// toes
//		0,
//		9.6730887889899994e-05,
//		3.0695446184800002e-12,
//		0,
//		{-4.1909515857699999e-09, 0, 0, 0},
//		0,
//		0
		};
//	double usrpos[3] = {35.16086834666666666, 139.6138257683333334, 47.347};
    double usrpos[3] = {35.48573997, 139.847435503, 1069679.0930026202};	
	double rs[6] = {0};
	double dopp = 0;
	double el=0.0;

	satEph_calc_satPos(time, &test_eph, rs);
	
	el = satellite_elevation_dopp(usrpos, rs, &dopp);
	
	
	TYP_F32 kk = SplitExtF64(dopp);
	
	printf("%.10f, %.10f\r\n", dopp - (double)kk.fDat - (double)kk.fErr, dopp - kk.fDat - kk.fErr);
	printf("dopp = %lf,  el = %lf\n", dopp, el);
    printf("Hello World \n");
	
	
	
	
//	el = satazel()
	
	
	
    return 0;
}

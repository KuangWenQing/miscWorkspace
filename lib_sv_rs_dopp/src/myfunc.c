#include "myfunc.h"

/* get satellite position , velocity and DOPP
 * 
 * args   : time_t time        I   second of week (s)
 *          eph_t   *eph       I   broadcast ephemeris
 *          double *usrpos     I   user position {lat,lon,high} (rad | m)
 *          double *rs         O   sat position and velocity (ecef)
 *                                 {x,y,z,vx,vy,vz} (m|m/s)
 *          double *dopp       O   satellite DOPP (Hz)
 * 
 * return : none
 * *************************************************/
void satellite_pos_vel_dopp(time_t time, const eph_t *eph, const double *usrpos, double *rs, double *dopp)
{
	STU_COOR_XYZ sat_vel, sat_pos, xyz;
	int i;
	double rst[3],tt=1E-3;
//	double dtst[1];

	eph2pos(time, eph, rs);
	time=timeadd(time,tt);
	eph2pos(time,eph,rst);
	
	
	/* satellite velocity by differential approx */
    for (i=0;i<3;i++) rs[i+3]=(rst[i]-rs[i])/tt;
	
	sat_vel.tX = rs[3];
	sat_vel.tY = rs[4];
	sat_vel.tZ = rs[5];
	
	sat_pos.tX = rs[0];
	sat_pos.tY = rs[1];
	sat_pos.tZ = rs[2];
	
	/* longitude latitude elevation  2  xyz */
	pos2ecef(usrpos, (double *) &xyz);
	/* calc DOPP */
	*dopp = calc_doppler_based_pv((PSTU_COOR_XYZ) &sat_vel, (PSTU_COOR_XYZ) &sat_pos, (PSTU_COOR_XYZ)&xyz);
}



double calc_doppler_based_pv(const PSTU_COOR_XYZ satVel, const PSTU_COOR_XYZ satPos, const PSTU_COOR_XYZ usrPos)
{
    double range, vs, cosU, tmpVall, tmpVal2;
    double deltaX, deltaY, deltaZ, tFqDivC;
    
    deltaX = satPos->tX - usrPos->tX;
    deltaY = satPos->tY - usrPos->tY;
    deltaZ = satPos->tZ - usrPos->tZ;
    
    // get the distance between sv and usr
    range = cal_range_3dim(satPos, usrPos);
    // get the velocity of satellite
    vs = cal_norm_3dim(satVel);     
    
    tmpVall = (deltaX * satVel->tX) + (deltaY * satVel->tY);
    tmpVal2 = deltaZ * satVel->tZ;
    tmpVall = tmpVall + tmpVal2;
    
    cosU = tmpVall / (range * vs); 
    tmpVall = (-1) * vs * cosU;
    
    tFqDivC = 5.25390625 + 0.00112915039;
    
    return tFqDivC * tmpVall;
}


double cal_range_3dim(const PSTU_COOR_XYZ Pos1, const PSTU_COOR_XYZ Pos2)
{
	return sqrt(pow((Pos1->tX - Pos2->tX), 2) + pow((Pos1->tY - Pos2->tY), 2) + pow((Pos1->tZ - Pos2->tZ), 2));
}

double cal_norm_3dim(const PSTU_COOR_XYZ V)
{
	return sqrt(pow(V->tX, 2) + pow(V->tY, 2) + pow(V->tZ, 2));
}

double timeadd(time_t t, double sec)
{
	time_t tmp;
	tmp = t + sec;
	
	if (tmp < -302400)
		tmp += 604800;
	else if(tmp > 302400)
		tmp -= 604800;
	
    return tmp;
}

void eph2pos(time_t time, const eph_t *eph, double *rs)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
    int n,sys,prn;
    
//    printf("eph2pos : time=%s sat=%2d\n",time_str(time,3),eph->sat);
    
    if (eph->A<=0.0) {
        rs[0]=rs[1]=rs[2]=0.0;
        return;
    }
    tk=timediff(time,eph->toe);
    
    switch ((sys=satsys(eph->sat,&prn))) {
        case SYS_GAL: mu=MU_GAL; omge=OMGE_GAL; break;
        case SYS_CMP: mu=MU_CMP; omge=OMGE_CMP; break;
        default:      mu=MU_GPS; omge=OMGE;     break;
    }
    M=eph->M0+(sqrt(mu/(eph->A*eph->A*eph->A))+eph->deln)*tk;
    
    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
        printf("kepler iteration overflow sat=%2d\n",eph->sat);
        return;
    }
    sinE=sin(E); cosE=cos(E);
    
    printf("kepler: sat=%2d e=%8.5f n=%2d del=%10.3e\n",eph->sat,eph->e,n,E-Ek);
    
    u=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;
    r=eph->A*(1.0-eph->e*cosE);
    i=eph->i0+eph->idot*tk;
    sin2u=sin(2.0*u); cos2u=cos(2.0*u);
    u+=eph->cus*sin2u+eph->cuc*cos2u;
    r+=eph->crs*sin2u+eph->crc*cos2u;
    i+=eph->cis*sin2u+eph->cic*cos2u;
    x=r*cos(u); y=r*sin(u); cosi=cos(i);
    
    /* beidou geo satellite (ref [9]) */
    if (sys==SYS_CMP&&prn<=5) {
        O=eph->OMG0+eph->OMGd*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        xg=x*cosO-y*cosi*sinO;
        yg=x*sinO+y*cosi*cosO;
        zg=y*sin(i);
        sino=sin(omge*tk); coso=cos(omge*tk);
        rs[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
        rs[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
        rs[2]=-yg*SIN_5+zg*COS_5;
    }
    else {
        O=eph->OMG0+(eph->OMGd-omge)*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        rs[0]=x*cosO-y*cosi*sinO;
        rs[1]=x*sinO+y*cosi*cosO;
        rs[2]=y*sin(i);
    }
}

time_t timediff(time_t t1, time_t t2)
{
//    return difftime(t1.time,t2.time)+t1.sec-t2.sec;
	time_t tmp;
	tmp = t1 - t2;
	
	if (tmp < -302400)
		tmp += 604800;
	else if(tmp > 302400)
		tmp -= 604800;
	
    return tmp;
}


int satsys(int sat, int *prn)
{
    int sys=SYS_NONE;
    if (sat<=0||MAXSAT<sat) sat=0;
    else if (sat<=NSATGPS) {
        sys=SYS_GPS; sat+=MINPRNGPS-1;
    }
    else if ((sat-=NSATGPS)<=NSATGLO) {
        sys=SYS_GLO; sat+=MINPRNGLO-1;
    }
    else if ((sat-=NSATGLO)<=NSATGAL) {
        sys=SYS_GAL; sat+=MINPRNGAL-1;
    }
    else if ((sat-=NSATGAL)<=NSATQZS) {
        sys=SYS_QZS; sat+=MINPRNQZS-1; 
    }
    else if ((sat-=NSATQZS)<=NSATCMP) {
        sys=SYS_CMP; sat+=MINPRNCMP-1; 
    }
    else if ((sat-=NSATCMP)<=NSATLEO) {
        sys=SYS_LEO; sat+=MINPRNLEO-1; 
    }
    else if ((sat-=NSATLEO)<=NSATSBS) {
        sys=SYS_SBS; sat+=MINPRNSBS-1; 
    }
    else sat=0;
    if (prn) *prn=sat;
    return sys;
}

/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos)
{
    double e2=FE_WGS84*(2.0-FE_WGS84),r2=dot(r,r,2),z,zk,v=RE_WGS84,sinp;
    
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}
/* transform geodetic to ecef position -----------------------------------------
* transform geodetic position to ecef position
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *r        O   ecef position {x,y,z} (m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void pos2ecef(const double *pos, double *r)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    double e2=FE_WGS84*(2.0-FE_WGS84),v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
    
    r[0]=(v+pos[2])*cosp*cosl;
    r[1]=(v+pos[2])*cosp*sinl;
    r[2]=(v*(1.0-e2)+pos[2])*sinp;
}


/* multiply matrix -----------------------------------------------------------*/
void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    double d;
    int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);
    
    for (i=0;i<n;i++) for (j=0;j<k;j++) {
        d=0.0;
        switch (f) {
            case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
            case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
            case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
            case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
        }
        if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
    }
}

void xyz2enu(const double *pos, double *E)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    
    E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
    E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
    E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
}

void ecef2enu(const double *pos, const double *r, double *e)
{
    double E[9];
    
    xyz2enu(pos,E);
    matmul("NN",3,1,3,1.0,E,r,0.0,e);
}

double dot(const double *a, const double *b, int n)
{
    double c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}

/* euclid norm -----------------------------------------------------------------
* euclid norm of vector
* args   : double *a        I   vector a (n x 1)
*          int    n         I   size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
extern double norm(const double *a, int n)
{
    return sqrt(dot(a,a,n));
}

/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position (ecef at transmission) (m)
*          double *rr       I   receiver position (ecef at reception) (m)
*          double *e        O   line-of-sight vector (ecef)
* return : geometric distance (m) (0>:error/no satellite position)
* notes  : distance includes sagnac effect correction
*-----------------------------------------------------------------------------*/
extern double geodist(const double *rs, const double *rr, double *e)
{
    double r;
    int i;
    
    if (norm(rs,3)<RE_WGS84) return -1.0;
    for (i=0;i<3;i++) e[i]=rs[i]-rr[i];
    r=norm(e,3);
    for (i=0;i<3;i++) e[i]/=r;
    return r+OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
}

/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite azimuth/elevation angle
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *e        I   receiver-to-satellilte unit vevtor (ecef)
*          double *azel     IO  azimuth/elevation {az,el} (rad) (NULL: no output)
*                               (0.0<=azel[0]<2*pi,-pi/2<=azel[1]<=pi/2)
* return : elevation angle (rad)
*-----------------------------------------------------------------------------*/
double satazel(const double *pos, const double *e, double *azel)
{
    double az=0.0,el=PI/2.0,enu[3];
    
    if (pos[2]>-RE_WGS84) {
        ecef2enu(pos,e,enu);
        az=dot(enu,enu,2)<1E-12?0.0:atan2(enu[0],enu[1]);
        if (az<0.0) az+=2*PI;
        el=asin(enu[2]);
    }
    if (azel) {azel[0]=az; azel[1]=el;}
    return el;
}


/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite elevation angle
* args   : double *usrpos   I   user position, geodetic position {lat,lon,h} (rad,m)
*          double *satpos   I   satellite position, 
* 
* return : elevation angle (rad)  -pi/2<= el <=pi/2  is right
* ----------------------------------------------*/
double satellite_elevation(const double *usrpos, const double* satpos)
{
	double e[3], rr[3], azel[2], r;
	
	pos2ecef(usrpos, rr);
	
	/* geometric distance/azimuth/elevation angle */
	if ((r=geodist(satpos, rr, e))<=0.0)
		return PI;
    return satazel(usrpos, e, azel);
}


TYP_F32 SplitExtF64(double dData)
{
	TYP_F32 tRes;
	double dTmp = dData * 0x20000001;
	
	tRes.fDat = dTmp - (dTmp - dData);
	tRes.fErr = dData - tRes.fDat;
	return tRes;
}

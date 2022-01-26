#include "myfunc.h"


/* get satellite position , velocity and DOPP
 * 
 * args   : gtime_t time       I   time (gpst)
 *          eph_t   *eph       I   broadcast ephemeris
 *          double *usrpos     I   user position
 *          double *rs         O   sat position and velocity (ecef)
 *                                 {x,y,z,vx,vy,vz} (m|m/s)
 *          double *dopp       O   satellite DOPP (Hz)
 * 
 * return : none
 * *************************************************/
void satellite_pos_vel_dopp(gtime_t time, const eph_t *eph, const double *usrpos, double *rs, double *dopp)
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

gtime_t timeadd(gtime_t t, double sec)
{
    double tt;
    
    t.sec+=sec; tt=floor(t.sec); t.time+=(int)tt; t.sec-=tt;
    return t;
}

void eph2pos(gtime_t time, const eph_t *eph, double *rs)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
    int n,sys,prn;
    
    printf("eph2pos : time=%s sat=%2d\n",time_str(time,3),eph->sat);
    
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

void time2epoch(gtime_t t, double *ep)
{
    const int mday[]={ /* # of days in a month */
        31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
        31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
    };
    int days,sec,mon,day;
    
    /* leap year if year%4==0 in 1901-2099 */
    days=(int)(t.time/86400);
    sec=(int)(t.time-(time_t)days*86400);
    for (day=days%1461,mon=0;mon<48;mon++) {
        if (day>=mday[mon]) day-=mday[mon]; else break;
    }
    ep[0]=1970+days/1461*4+mon/12; ep[1]=mon%12+1; ep[2]=day+1;
    ep[3]=sec/3600; ep[4]=sec%3600/60; ep[5]=sec%60+t.sec;
}

void time2str(gtime_t t, char *s, int n)
{
    double ep[6];
    
    if (n<0) n=0; else if (n>12) n=12;
    if (1.0-t.sec<0.5/pow(10.0,n)) {t.time++; t.sec=0.0;};
    time2epoch(t,ep);
    sprintf(s,"%04.0f/%02.0f/%02.0f %02.0f:%02.0f:%0*.*f",ep[0],ep[1],ep[2],
            ep[3],ep[4],n<=0?2:n+3,n<=0?0:n,ep[5]);
}

char *time_str(gtime_t t, int n)
{
    static char buff[64];
    time2str(t,buff,n);
    return buff;
}


double timediff(gtime_t t1, gtime_t t2)
{
    return difftime(t1.time,t2.time)+t1.sec-t2.sec;
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



/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite azimuth/elevation angle
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *e        I   receiver-to-satellilte unit vevtor (ecef)
*          double *azel     IO  azimuth/elevation {az,el} (rad) (NULL: no output)
*                               (0.0<=azel[0]<2*pi,-pi/2<=azel[1]<=pi/2)
* return : elevation angle (rad)
*-----------------------------------------------------------------------------*/
//double satazel(const double *pos, const double *e, double *azel)
//{
//    double az=0.0,el=PI/2.0,enu[3];
//    
//    if (pos[2]>-RE_WGS84) {
//        ecef2enu(pos,e,enu);
//        az=dot(enu,enu,2)<1E-12?0.0:atan2(enu[0],enu[1]);
//        if (az<0.0) az+=2*PI;
//        el=asin(enu[2]);
//    }
//    if (azel) {azel[0]=az; azel[1]=el;}
//    return el;
//}
//
//
//void matmul(const char *tr, int n, int k, int m, double alpha,
//                   const double *A, const double *B, double beta, double *C)
//{
//    int lda=tr[0]=='T'?m:n,ldb=tr[1]=='T'?k:m;
//    
//    dgemm_((char *)tr,(char *)tr+1,&n,&k,&m,&alpha,(double *)A,&lda,(double *)B,
//           &ldb,&beta,C,&n);
//}
//
//void xyz2enu(const double *pos, double *E)
//{
//    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
//    
//    E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
//    E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
//    E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
//}
//
//void ecef2enu(const double *pos, const double *r, double *e)
//{
//    double E[9];
//    
//    xyz2enu(pos,E);
//    matmul("NN",3,1,3,1.0,E,r,0.0,e);
//}
//
//
double dot(const double *a, const double *b, int n)
{
    double c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}
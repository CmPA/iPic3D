/*******************************************************************************************
EllipticF.cpp  -  Elliptic Functions Tool
                            -------------------
developers: Found on the web
********************************************************************************************/
#include "EllipticF.h"

#define MAX(x,y) (x>y?x:y)
#define MAX3(x,y,z) MAX(MAX(x,y),z)
#define MIN(x,y) (x>y?y:x)
#define MIN3(x,y,z) MIN(MIN(x,y),z)


/*
 * NOTE: constants declared in float.h for the x86 CPU are:
 *
 * DBL_EPSILON	2.2204460492503131e-016,smallest such that 1.0+DBL_EPSILON!=1.0
 * DBL_MAX 	1.7976931348623158e+308  maximum possible value for type double
 * DBL_MIN 	2.2250738585072014e-308  minimum possible positive value for double
 *
 * If you are compiling this code for a non-Intel CPU, your values for these
 * constants may be different.
 */



/*
 * EllipticF(k): complete elliptic integral of the first kind
 */

//#define EllipticF(k,ierr) drf(0.0,1.0-pow(k,2),1.0,&ierr)

/*
 * drf.c:   Compute the complete or incomplete elliptic integral of the
 *          first kind.
 *
 * Description:
 *
 *  For x,y and z non-negative and at most one of them zero, drf(x,y,z)
 *  = integral from zero to infinity of 
 *
 *
 *            -1/2     -1/2     -1/2
 *  (1/2)(t+x)    (t+y)    (t+z)    dt.
 *
 *  If x, y or z is zero, the integral is complete.
 *
 *  *piErr returns non-zero if any arguments are invalid.
 *
 * Credits:
 *
 *  This function is adapted by E. Dennison from FORTRAN code written by:
 *
 *  Carlson, B.C.
 *  Notis, E.M
 *  Pexton, R.L.
 *
 */

double drf(double x, double y, double z, int* piErr)
{
int    iErr=0;
double mu,
       xn,yn,zn,
       xndev,yndev,zndev,
       xnroot,ynroot,znroot,
       lambda,
       epsilon,
       e2,e3,
       result,
       s;
	
    const double c1 = 1.0/24.0;
    const double c2 = 3.0/44.0;
    const double c3 = 1.0/14.0;
    const double errtol = pow(DBL_EPSILON*4.0,1.0/6.0);
    const double lolim = 5.0*DBL_MIN;
    const double hilim = DBL_MAX/5.0;

    if (piErr)
    {
        if (MIN3(x,y,z)<0.0)
        {
            iErr = 1;
        }
        else if (MIN3(x+y,x+z,y+z)<lolim)
        {
            iErr = 2;
        }
        else if (MAX3(x,y,z)>hilim)
        {
            iErr = 3;
        }
    }
    if (iErr)
    {
        if (piErr)
        {
            *piErr = iErr;
        }
        result = 0.0;
    }

else
    {
        xn = x;
        yn = y;
        zn = z;

        while (1)
        {
            mu = (xn+yn+zn)/3.0;
            xndev = 2.0-(mu+xn)/mu;
            yndev = 2.0-(mu+yn)/mu;
            zndev = 2.0-(mu+zn)/mu;
            epsilon = MAX3(fabs(xndev),fabs(yndev),fabs(zndev));
            if (epsilon < errtol) break;
            xnroot = sqrt(xn);
            ynroot = sqrt(yn);
            znroot = sqrt(zn);
            lambda = xnroot*(ynroot+znroot) +ynroot*znroot;
            xn = (xn+lambda)*0.25;
            yn = (yn+lambda)*0.25;
            zn = (zn+lambda)*0.25;
        }
        e2 = xndev*yndev - pow(zndev,2);
        e3 = xndev*yndev*zndev;
        s = 1.0 + (c1*e2-0.1-c2*e3)*e2 + c3*e3;

        if (piErr)
        {
            *piErr = 0;
        }
        result = s/sqrt(mu);
    }
    return result;
}
/* END drf() */


/*
 * E(k): complete elliptic integral of the second kind
 */

//#define EllipticE(k,ierr)   (drf(0.0,1.0-pow(k,2),1.0,&ierr)\
                    -(pow(k,2)/3.0)*drd(0.0,1.0-pow(k,2),1.0,&ierr))

 /*
  * FastE(k): fast, complete elliptic integral of the second kind.
  *           Use this macro if the complete elliptic integral of
  *           the first kind was previously computed for the same
  *           value of k.
  */

#define FastE(F,k,ierr)	((EllipticF)-(pow(k,2)/3.0)*drd(0.0,1.0-pow(k,2),1.0,&ierr))

/*
 * drd.c:   Compute the complete or incomplete elliptic integral of the
 *          second kind.
 *
 * Description:
 *
 *  For x and y non-negative, x+y and z positive, drf(x,y,z) = integral 
 *  from zero to infinity of 
 *
 *
 *            -1/2     -1/2     -3/2
 *  (3/2)(t+x)    (t+y)    (t+z)    dt.
 *
 *  If x or y is zero, the integral is complete.
 *
 *  *piErr returns non-zero if any arguments are invalid.
 *
 * Credits:
 *
 *  This function is adapted by E. Dennison from FORTRAN code written by:
 *
 *  Carlson, B.C.
 *  Notis, E.M
 *  Pexton, R.L.
 *
 */

double drd(double x, double y, double z, int* piErr)
{
    int     iErr=0;
    double  mu,
            xn,yn,zn,
            xndev,yndev,zndev,
            xnroot,ynroot,znroot,
            lambda,
            epsilon,
            ea,eb,ec,ed,ef,
            sigma,
            power4,
            result,
            s1,s2;
	
    const double c1 = 3.0/14.0;
    const double c2 = 1.0/6.0;
    const double c3 = 9.0/22.0;
    const double c4 = 3.0/26.0;
    const double errtol = pow(DBL_EPSILON/3.0,1.0/6.0);
    double uplim;
    const double lolim = 2.0/pow(DBL_MAX,2.0/3.0);
    double tuplim = pow(DBL_MIN,1.0/3.0);
    tuplim = pow(0.1*errtol,1.0/3.0)/tuplim;
    uplim = pow(tuplim,2.0);

    if (piErr)
    {
        if (MIN(x,y)<0.0)
        {
            iErr = 1;
        }
        else if (MAX3(x,y,z)>uplim)
        {
            iErr = 2;
        }
        else if (MIN(x+y,z)<lolim)
        {
            iErr = 3;
        }
    }
    if (iErr)
    {
        if (piErr)
        {
            *piErr = iErr;
        }
        result = 0.0;
    }

else
    {
        xn = x;
        yn = y;
        zn = z;
        sigma = 0.0;
        power4 = 1.0;
        while (1)
        {
            mu = (xn+yn+3.0*zn)*0.2;
            xndev = (mu-xn)/mu;
            yndev = (mu-yn)/mu;
            zndev = (mu-zn)/mu;
            epsilon = MAX3(fabs(xndev),fabs(yndev),fabs(zndev));
            if (epsilon < errtol) break;
            xnroot = sqrt(xn);
            ynroot = sqrt(yn);
            znroot = sqrt(zn);
            lambda = xnroot*(ynroot+znroot) +ynroot*znroot;
            sigma = sigma+power4/(znroot*(zn+lambda));
            power4 = power4*0.25;
            xn = (xn+lambda)*0.25;
            yn = (yn+lambda)*0.25;
            zn = (zn+lambda)*0.25;
        }
        ea = xndev*yndev;
        eb = zndev*zndev;
        ec = ea-eb;
        ed = ea-6.0*eb;
        ef = ed+ec+ec;
        s1 = ed*(-c1+0.25*c3*ed-1.5*c4*zndev*ef);
        s2 = zndev*(c2*ef+zndev*(-c3*ec+zndev*c4*ea));
        if (piErr)
        {
            *piErr = 0;
        }
        result = 3.0*sigma+power4*(1.0+s1+s2)/(mu*sqrt(mu));
    }
    return result;
}

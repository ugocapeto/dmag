#include "header.h"
#include "proto.h"

void psi_prime_gradient(
 double *I1,
 double *I1x,
 double *I1y,
 double *I2w,
 double *I2wx,
 double *I2wy,
 double *I2wxx,
 double *I2wxy,
 double *I2wyy,
 double *du,
 double *dv,
 double *psig,
 int xdim,
 int ydim,
 double epsilon
)

/*
psig must be allocated outside this function
*/

{

 int i;
 int j;
 double val;
 double valx;
 double valy;
 double s2;

 for ( i= 0 ; i< ydim ; i++ ) { 
    for ( j= 0 ; j< xdim ; j++ ) {
       valx=
        I2wx[i*xdim+j]+
        I2wxx[i*xdim+j]*du[i*xdim+j]+
        I2wxy[i*xdim+j]*dv[i*xdim+j]-
        I1x[i*xdim+j];
       valy=
        I2wy[i*xdim+j]+
        I2wxy[i*xdim+j]*du[i*xdim+j]+
        I2wyy[i*xdim+j]*dv[i*xdim+j]-
        I1y[i*xdim+j];
       s2= valx*valx+valy*valy;
       val= sqrt(s2+epsilon*epsilon);
       val= 1.0/(2.0*val);
       psig[i*xdim+j]= val;
    }
 }

}

#include "header.h"
#include "proto.h"

void psi_prime_data(
 double *I1,
 double *I2w,
 double *I2wx,
 double *I2wy,
 double *du,
 double *dv,
 double *psid,
 int xdim,
 int ydim,
 double epsilon
)

/*
psid must be allocated outside this function
*/

{

 int i;
 int j;
 double val;
 double s2;

 for ( i= 0 ; i< ydim ; i++ ) { 
    for ( j= 0 ; j< xdim ; j++ ) {
       val=
        I2w[i*xdim+j]+
        I2wx[i*xdim+j]*du[i*xdim+j]+
        I2wy[i*xdim+j]*dv[i*xdim+j]-
        I1[i*xdim+j];
       s2= val*val;
       val= sqrt(s2+epsilon*epsilon);
       val= 1.0/(2.0*val);
       psid[i*xdim+j]= val;
    }
 }

}

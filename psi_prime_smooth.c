#include "header.h"
#include "proto.h"

void psi_prime_smooth(
 double *ux,
 double *uy,
 double *vx,
 double *vy,
 double *psis,
 int xdim,
 int ydim,
 double epsilon
)

/*
psis must be allocated outside this function
*/

{

 int i;
 int j;
 double val;
 double valu;
 double valv;
 double s2;

 for ( i= 0 ; i< ydim ; i++ ) { 
    for ( j= 0 ; j< xdim ; j++ ) {
       valu= ux[i*xdim+j]*ux[i*xdim+j]+
             uy[i*xdim+j]*uy[i*xdim+j];
       valv= vx[i*xdim+j]*vx[i*xdim+j]+
             vy[i*xdim+j]*vy[i*xdim+j];
       s2= valu+valv; 
       val= sqrt(s2+epsilon*epsilon);
       val= 1.0/(2.0*val);
       psis[i*xdim+j]= val;
    }
 }

}

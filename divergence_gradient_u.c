#include "header.h"
#include "proto.h"

void divergence_gradient_u(
 double *u,
 double *psi1,
 double *psi2,
 double *psi3,
 double *psi4,
 double *div_u,
 int xdim,
 int ydim
)

/*
div_u must be allocated outside this function
*/

/*
In the paper,
i is associated with x and j is associated with y
but in the code,
i is associated with y and j is associated with x
which means that
i must be swapped with j and j must be swapped with i
*/

{

 int i;
 int j;

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {

       div_u[i*xdim+j]= 0.0;

       if ( (j+1) < xdim ) {
          div_u[i*xdim+j]+=
           psi1[i*xdim+j]*(u[i*xdim+(j+1)]-u[i*xdim+j]);
       }

       if ( (j-1) >= 0 ) {
          div_u[i*xdim+j]+=
           psi2[i*xdim+j]*(u[i*xdim+(j-1)]-u[i*xdim+j]);
       }

       if ( (i+1) < ydim ) {
          div_u[i*xdim+j]+=
           psi3[i*xdim+j]*(u[(i+1)*xdim+j]-u[i*xdim+j]);
       }

       if ( (i-1) >= 0 ) {
          div_u[i*xdim+j]+=
           psi4[i*xdim+j]*(u[(i-1)*xdim+j]-u[i*xdim+j]);
       }
    }
 }

}

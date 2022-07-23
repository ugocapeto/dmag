#include "header.h"
#include "proto.h"

void divergence_gradient_v(
 double *v,
 double *psi1,
 double *psi2,
 double *psi3,
 double *psi4,
 double *div_v,
 int xdim,
 int ydim
)

/*
div_v must be allocated outside this function
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

       div_v[i*xdim+j]= 0.0;

       if ( (j+1) < xdim ) {
          div_v[i*xdim+j]+=
           psi1[i*xdim+j]*(v[i*xdim+(j+1)]-v[i*xdim+j]);
       }
       if ( (j-1) >= 0 ) {
          div_v[i*xdim+j]+=
           psi2[i*xdim+j]*(v[i*xdim+(j-1)]-v[i*xdim+j]);
       }

       if ( (i+1) < ydim ) {
          div_v[i*xdim+j]+=
           psi3[i*xdim+j]*(v[(i+1)*xdim+j]-v[i*xdim+j]);
       }

       if ( (i-1) >= 0 ) {
          div_v[i*xdim+j]+=
           psi4[i*xdim+j]*(v[(i-1)*xdim+j]-v[i*xdim+j]);
       }
    }
 }

}

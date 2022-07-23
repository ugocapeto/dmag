#include "header.h"
#include "proto.h"

void psi_prime_smooth_1234(
 double *psis,
 double *psi1,
 double *psi2,
 double *psi3,
 double *psi4,
 int xdim,
 int ydim
)

/*
psi1,psi2,psi3,psi4 must be allocated outside this function
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

       if ( (j+1) < xdim )
        psi1[i*xdim+j]= 0.5*(psis[i*xdim+(j+1)]+psis[i*xdim+j]);
       else
        psi1[i*xdim+j]= psis[i*xdim+j];

       if ( (j-1) >= 0 )
        psi2[i*xdim+j]= 0.5*(psis[i*xdim+(j-1)]+psis[i*xdim+j]);
       else
        psi2[i*xdim+j]= psis[i*xdim+j];

       if ( (i+1) < ydim )
        psi3[i*xdim+j]= 0.5*(psis[(i+1)*xdim+j]+psis[i*xdim+j]);
       else
        psi3[i*xdim+j]= psis[i*xdim+j];

       if ( (i-1) >= 0 )
        psi4[i*xdim+j]= 0.5*(psis[(i-1)*xdim+j]+psis[i*xdim+j]);
       else
        psi4[i*xdim+j]= psis[i*xdim+j];

    }
 }

}

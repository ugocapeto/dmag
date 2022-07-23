#include "header.h"
#include "proto.h"

void image_derivative_x(
 double *I,
 double *Ix,
 int xdim,
 int ydim
)

/*
Compute the image derivative in the x direction
*/

/*
The image derivative should be allocated outside
*/

{

 int i;
 int j;

 /*
 Compute the derivative using central differences
 except at column 0 and column (xdim-1)
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 1 ; j< (xdim-1) ; j++ ) {
       Ix[xdim*i+j]= 0.5*(I[xdim*i+(j+1)]-I[xdim*i+(j-1)]);
    }
 }

 /*
 Compute the derivative using forward differences
 at column 0
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    j= 0;
    Ix[xdim*i+j]= I[xdim*i+(j+1)]-I[xdim*i+j];
 }

 /*
 Compute the derivative using backward differences
 at column (xdim-1)
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    j= xdim-1;
    Ix[xdim*i+j]= I[xdim*i+j]-I[xdim*i+(j-1)];
 }

}

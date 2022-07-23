#include "header.h"
#include "proto.h"

void image_derivative_y(
 double *I,
 double *Iy,
 int xdim,
 int ydim
)

/*
Compute the image derivative in the y direction
*/

/*
The image derivative should be allocated outside
*/

{

 int i;
 int j;

 /*
 Compute the derivative using central differences
 except at row 0 and row (ydim-1)
 */

 for ( i= 1 ; i< (ydim-1) ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       Iy[xdim*i+j]= 0.5*(I[xdim*(i+1)+j]-I[xdim*(i-1)+j]);
    }
 }

 /*
 Compute the derivative using forward differences
 at row 0
 */

 i= 0;
 for ( j= 0 ; j< xdim ; j++ ) {
    Iy[xdim*i+j]= I[xdim*(i+1)+j]-I[xdim*i+j];
 }

 /*
 Compute the derivative using backward differences
 at row (ydim-1)
 */

 i= ydim-1;
 for ( j= 0 ; j< xdim ; j++ ) {
    Iy[xdim*i+j]= I[xdim*i+j]-I[xdim*(i-1)+j];
 }

}

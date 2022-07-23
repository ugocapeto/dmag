#include "header.h"
#include "proto.h"

void warp_image(
 double *inp_I,
 double *u,
 double *v,
 double *out_I,
 int xdim,
 int ydim
)

/*
Warp the input image using the u and v displacement arrays
*/

/*
The output image should be allocated outside
*/

{

 int i;
 int j;
 double disp_u;
 double disp_v;
 double x;
 double y;
 double f;

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {

       /*
       Get the displacement vector (w) at that pixel
       */

       disp_u= u[xdim*i+j];
       disp_v= v[xdim*i+j];

       /* 
       Get the offset location we need to consider
       */

       x= (double)j+disp_u;
       y= (double)i+disp_v;

       /*
       Make sure x is in bounds before interpolating
       */

       /*
       if ( x < 0.0 )
        x= 0.0;
       if ( x > (double)(xdim-1) )
        x= (double)(xdim-1);
       */
       if ( x < 0.0 ) {
          out_I[xdim*i+j]= 0.0;
          continue;
       }
       if ( x > (double)(xdim-1) ) {
          out_I[xdim*i+j]= 0.0;
          continue;
       }

       /*
       Make sure y is in bounds before interpolating
       */

       /* 
       if ( y < 0.0 )
        y= 0.0;
       if ( y > (double)(ydim-1) )
        y= (double)(ydim-1);
       */
       if ( y < 0.0 ) {
          out_I[xdim*i+j]= 0.0;
          continue;
       }
       if ( y > (double)(ydim-1) ) {
          out_I[xdim*i+j]= 0.0;
          continue;
       }

       /*
       Interpolate 
       */

       f= bicubic_interpolation_on_image(
        inp_I,
        x,
        y,
        xdim,
        ydim
       );

       out_I[xdim*i+j]= f; 
    }
 }

}

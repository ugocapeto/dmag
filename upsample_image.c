#include "header.h"
#include "proto.h"

void upsample_image(
 double *inp_I,
 int inp_xdim,
 int inp_ydim,
 double *out_I,
 int out_xdim,
 int out_ydim
)

/*
The output image must be allocated outside the function
*/

{

 double factor_x;
 double factor_y;
 int i;
 int j;
 double x;
 double y;
 double f;

 /*
 Compute the zoom factor in both directions
 */

 factor_x= (double)out_xdim/(double)inp_xdim;
 factor_y= (double)out_ydim/(double)inp_ydim;

 for ( i= 0 ; i< out_ydim ; i++ ) {
    for ( j= 0 ; j< out_xdim ; j++ ) {

       /*
       Get the corresponding pixel location
       in the input image
       */

       x= (double)j/factor_x;
       y= (double)i/factor_y;

       /*
       Make sure (x,y) are in bounds
       */

       if ( x < 0.0 )
        x= 0.0;
       if ( x > (inp_xdim-1) )
        x= (inp_xdim-1);
       if ( y < 0.0 )
        y= 0.0;
       if ( y > (inp_ydim-1) )
        y= (inp_ydim-1);

       /*
       Interpolate the value using bicubic interpolation
       */

       f= bicubic_interpolation_on_image(
        inp_I,
        x,
        y,
        inp_xdim,
        inp_ydim
       );

       out_I[out_xdim*i+j]= f;
    }
 }

}

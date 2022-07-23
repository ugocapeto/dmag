#include "header.h"
#include "proto.h"

void downsample_image(
 double *inp_I,
 int inp_xdim,
 int inp_ydim,
 double *out_I,
 int out_xdim,
 int out_ydim,
 double factor
)

/*
The output image must be allocated outside the function
*/

{

 int size;
 double *I;
 int ind;
 double zoom_sigma_zero;
 double sigma;
 int precision;
 int i;
 int j;
 double xval;
 double yval;
 double f;

 /*
 Copy the input image
 */

 size= inp_xdim*inp_ydim;
 I= (double *)calloc(size,sizeof(double));
 for ( ind= 0 ; ind< size ; ind++ )
  I[ind]= inp_I[ind];

 /*
 Apply a Gaussian blur to the copied image
 */

 zoom_sigma_zero= 0.6;
 sigma= zoom_sigma_zero*sqrt(1./(factor*factor)-1.);

 precision= 5;
 gaussian_blur_image(
  I,
  inp_xdim,
  inp_ydim,
  sigma,
  precision
 );

 for ( i= 0 ; i< out_ydim ; i++ ) {
    for ( j= 0 ; j< out_xdim ; j++ ) {

       /*
       Get the coordinates in the blurred image
       */

       yval= (double)i/factor;
       xval= (double)j/factor;

       if ( yval < 0 )
        yval= 0;
       if ( yval > (inp_ydim-1) )
        yval= inp_ydim-1;
       if ( xval < 0 )
        xval= 0;
       if ( xval > (inp_xdim-1) )
        xval= inp_xdim-1;

       /*
       Interpolate the intensity
       */

       f= bicubic_interpolation_on_image(
        I,
        xval,
        yval,
        inp_xdim,
        inp_ydim
       );

       out_I[out_xdim*i+j]= f;
    }
 }

 /*
 Free the copy of the input image
 */

 free(I);

}

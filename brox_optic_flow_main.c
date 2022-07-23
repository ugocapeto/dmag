#include "header.h"
#include "proto.h"

void brox_optic_flow_main(
 double *I1,
 double *I2,
 double *disp_u,
 double *disp_v,
 int w,
 int h,
 double factor,
 int nscales,
 double epsilon,
 double gamma,
 double alpha,
 int max_inner_iter,
 int max_outer_iter,
 double sor_w,
 int max_sor_iter,
 double max_sor_error
)

/*
disp_u,disp_v are output displacements
*/

{

 int nx;
 int ny;
 double **I1s;
 double **I2s;
 int *nxs;
 int *nys;
 double **us;
 double **vs;
 int size;
 double sigma;
 int precision;
 int s;
 int i;
 int j;
 int xdim;
 int ydim;

 /*
 We are gonna store the images and displacements
 at each scale
 0           corresponds to the finest scale (bottom of pyramid)
 (nscales-1) corresponds to the coarsest scale (top of pyramid)
 */

 nx= w;
 ny= h;

 I1s= (double **)calloc(nscales,sizeof(double *));
 I2s= (double **)calloc(nscales,sizeof(double *));
 nxs= (int *)calloc(nscales,sizeof(int));
 nys= (int *)calloc(nscales,sizeof(int));
 us= (double **)calloc(nscales,sizeof(double *));
 vs= (double **)calloc(nscales,sizeof(double *));

 /*
 Let's take care of scale 0 (finest scale)
 */

 nxs[0]= nx;
 nys[0]= ny;
 size= nxs[0]*nys[0];
 I1s[0]= (double *)calloc(size,sizeof(double));
 I2s[0]= (double *)calloc(size,sizeof(double));
 us[0]= (double *)calloc(size,sizeof(double));
 vs[0]= (double *)calloc(size,sizeof(double));

 /*
 Let's normalize the input images
 and put the normalized images in I1s[0],I2s[0]
 */

 normalize_image_pair(
  I1,
  I2,
  I1s[0],
  I2s[0],
  size
 );

 /*
 Apply a Gaussian blur to image 1 at the finest scale
 */

 sigma= 0.8;
 precision= 5;
 gaussian_blur_image(
  I1s[0],
  nxs[0],
  nys[0],
  sigma,
  precision
 );

 /*
 Apply a Gaussian blur to image 2 at the finest scale
 */

 sigma= 0.8;
 precision= 5;
 gaussian_blur_image(
  I2s[0],
  nxs[0],
  nys[0],
  sigma,
  precision
 );

 /*
 Copy the input displacement maps into scale 0
 */

 xdim= nxs[0];
 ydim= nys[0];

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       us[0][i*xdim+j]= disp_u[i*xdim+j];
       vs[0][i*xdim+j]= disp_v[i*xdim+j];
    }
 }

 /*
 Let's build the pyramid of images
 by downsampling the input images
 */

 for ( s= 1 ; s< nscales ; s++ ) {

    /*
    Let's get the image dimensions at that scale
    using the dimensions from the previous scale
    */

    downsample_image_dimensions(
     nxs[s-1],
     nys[s-1],
     &nxs[s],
     &nys[s],
     factor
    );

    size= nxs[s]*nys[s];
    I1s[s]= (double *)calloc(size,sizeof(double));
    I2s[s]= (double *)calloc(size,sizeof(double));
    us[s]= (double *)calloc(size,sizeof(double));
    vs[s]= (double *)calloc(size,sizeof(double));

    /*
    Downsample image 1
    */

    downsample_image(
     I1s[s-1],
     nxs[s-1],
     nys[s-1],
     I1s[s],
     nxs[s],
     nys[s],
     factor
    );

    /*
    Downsample image 2
    */

    downsample_image(
     I2s[s-1],
     nxs[s-1],
     nys[s-1],
     I2s[s],
     nxs[s],
     nys[s],
     factor
    );

 }

 /*
 It's time to actually do something
 */

 for ( s= nscales-1 ; s>= 0 ; s-- ) {

    /*
    Let's update the displacement at scale s
    */

    fprintf(stdout,"Processing at scale = %d (width = %d height = %d)\n",s,nxs[s],nys[s]);
	
    brox_optic_flow(
     I1s[s],
     I2s[s],
     us[s],
     vs[s],
     nxs[s],
     nys[s],
     epsilon,
     gamma,
     alpha,
     max_inner_iter,
     max_outer_iter,
     sor_w,
     max_sor_iter,
     max_sor_error
    );

    if ( s == 0 ) {

       /*
       We are at the finer scale, so we are done
       */

       break;
    }

    /*
    Upsample the displacement for the next finer scale
    */

    upsample_image(
     us[s],
     nxs[s],
     nys[s],
     us[s-1],
     nxs[s-1],
     nys[s-1]
    );

    upsample_image(
     vs[s],
     nxs[s],
     nys[s],
     vs[s-1],
     nxs[s-1],
     nys[s-1]
    );

    /*
    Displacement values
    should be divided by the downsampling factor
    */

    for ( i= 0 ; i< nys[s-1] ; i++ ) {
       for ( j= 0 ; j< nxs[s-1] ; j++ ) {
          us[s-1][nxs[s-1]*i+j]/= factor;
          vs[s-1][nxs[s-1]*i+j]/= factor;
       }
    }
 }

 /*
 Dump the displacement maps
 */

 xdim= nxs[0];
 ydim= nys[0];

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       disp_u[i*xdim+j]= us[0][i*xdim+j];
       disp_v[i*xdim+j]= vs[0][i*xdim+j];
    }
 }

 /*
 Free the pyramid of images and displacements
 */

 for ( s= 0 ; s< nscales ; s++ ) {
    free(I1s[s]);
    free(I2s[s]);
    free(us[s]);
    free(vs[s]);
 }
 free(I1s);
 free(I2s);
 free(nxs);
 free(nys);
 free(us);
 free(vs);
 
}

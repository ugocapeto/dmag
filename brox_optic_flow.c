#include "header.h"
#include "proto.h"

void brox_optic_flow(
 double *I1,
 double *I2,
 double *u,
 double *v,
 int xdim,
 int ydim,
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
The displacements u,v are both input and output
*/

{

 double *I1x;
 double *I1y;
 double *I2x;
 double *I2y;
 double *I2xx;
 double *I2xy;
 double *I2yy;

 /*
 Allocate memory
 */

 I1x= (double *)calloc(xdim*ydim,sizeof(double));
 I1y= (double *)calloc(xdim*ydim,sizeof(double));
 I2x= (double *)calloc(xdim*ydim,sizeof(double));
 I2y= (double *)calloc(xdim*ydim,sizeof(double));
 I2xx= (double *)calloc(xdim*ydim,sizeof(double));
 I2xy= (double *)calloc(xdim*ydim,sizeof(double));
 I2yy= (double *)calloc(xdim*ydim,sizeof(double));

 /*
 Compute first order derivatives for image 1
 */

 image_derivative_x(
  I1,
  I1x,
  xdim,
  ydim
 );

 image_derivative_y(
  I1,
  I1y,
  xdim,
  ydim
 );

 /*
 Compute first order derivatives for image 2
 */

 image_derivative_x(
  I2,
  I2x,
  xdim,
  ydim
 );

 image_derivative_y(
  I2,
  I2y,
  xdim,
  ydim
 );

 /*
 Compute second order derivatives for image 2
 */

 image_derivative_x(
  I2x,
  I2xx,
  xdim,
  ydim
 );

 image_derivative_y(
  I2x,
  I2xy,
  xdim,
  ydim
 );

 image_derivative_y(
  I2y,
  I2yy,
  xdim,
  ydim
 );

 /*
 Update the displacement
 */

 update_displacement(
  I1,
  I2,
  I1x,
  I1y,
  I2x,
  I2y,
  I2xx,
  I2xy,
  I2yy,
  u,
  v,
  xdim,
  ydim,
  epsilon,
  gamma,
  alpha,
  max_inner_iter,
  max_outer_iter,
  sor_w,
  max_sor_iter,
  max_sor_error
 );

 /*
 Free all the allocated memory
 */

 free(I1x);
 free(I1y);
 free(I2x);
 free(I2y);
 free(I2xx);
 free(I2xy);
 free(I2yy);

}

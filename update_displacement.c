#include "header.h"
#include "proto.h"

void update_displacement(
 double *I1,
 double *I2,
 double *I1x,
 double *I1y,
 double *I2x,
 double *I2y,
 double *I2xx,
 double *I2xy,
 double *I2yy,
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

 double *I2w;
 double *I2wx;
 double *I2wy;
 double *I2wxx;
 double *I2wxy;
 double *I2wyy;
 double *ux;
 double *uy;
 double *vx;
 double *vy;
 double *psis;
 double *psi1;
 double *psi2;
 double *psi3;
 double *psi4;
 double *div_u;
 double *div_v;
 double *div_d;
 double *du;
 double *dv;
 double *psid;
 double *psig;
 double *Au;
 double *Av;
 double *Du;
 double *Dv;
 double *D;
 int outer_iter;
 int i;
 int j;
 int inner_iter;
 double min_du_val;
 double max_du_val;
 double min_dv_val;
 double max_dv_val;
 double avg_du_val;
 double avg_dv_val;
 double du_val;
 double dv_val;

 /*
 Allocate memory
 */

 I2w= (double *)calloc(xdim*ydim,sizeof(double));
 I2wx= (double *)calloc(xdim*ydim,sizeof(double));
 I2wy= (double *)calloc(xdim*ydim,sizeof(double));
 I2wxx= (double *)calloc(xdim*ydim,sizeof(double));
 I2wxy= (double *)calloc(xdim*ydim,sizeof(double));
 I2wyy= (double *)calloc(xdim*ydim,sizeof(double));
 ux= (double *)calloc(xdim*ydim,sizeof(double));
 uy= (double *)calloc(xdim*ydim,sizeof(double));
 vx= (double *)calloc(xdim*ydim,sizeof(double));
 vy= (double *)calloc(xdim*ydim,sizeof(double));
 psis= (double *)calloc(xdim*ydim,sizeof(double));
 psi1= (double *)calloc(xdim*ydim,sizeof(double));
 psi2= (double *)calloc(xdim*ydim,sizeof(double));
 psi3= (double *)calloc(xdim*ydim,sizeof(double));
 psi4= (double *)calloc(xdim*ydim,sizeof(double));
 div_u= (double *)calloc(xdim*ydim,sizeof(double));
 div_v= (double *)calloc(xdim*ydim,sizeof(double));
 div_d= (double *)calloc(xdim*ydim,sizeof(double));
 du= (double *)calloc(xdim*ydim,sizeof(double));
 dv= (double *)calloc(xdim*ydim,sizeof(double));
 psid= (double *)calloc(xdim*ydim,sizeof(double));
 psig= (double *)calloc(xdim*ydim,sizeof(double));
 Au= (double *)calloc(xdim*ydim,sizeof(double));
 Av= (double *)calloc(xdim*ydim,sizeof(double));
 Du= (double *)calloc(xdim*ydim,sizeof(double));
 Dv= (double *)calloc(xdim*ydim,sizeof(double));
 D= (double *)calloc(xdim*ydim,sizeof(double));

 outer_iter= 0;

 START:

 /*
 Outer iteration
 */

 outer_iter++;
 
 fprintf(stdout,"Iteration (outer) = %d\n",outer_iter);

 /*
 Warp image 2 and its first and second order derivatives
 */

 warp_image(
  I2,
  u,
  v,
  I2w,
  xdim,
  ydim
 );

 warp_image(
  I2x,
  u,
  v,
  I2wx,
  xdim,
  ydim
 );

 warp_image(
  I2y,
  u,
  v,
  I2wy,
  xdim,
  ydim
 );

 warp_image(
  I2xx,
  u,
  v,
  I2wxx,
  xdim,
  ydim
 );

 warp_image(
  I2xy,
  u,
  v,
  I2wxy,
  xdim,
  ydim
 );

 warp_image(
  I2yy,
  u,
  v,
  I2wyy,
  xdim,
  ydim
 );

 /*
 Compute first order derivatives of the displacement 
 */

 image_derivative_x(
  u,
  ux,
  xdim,
  ydim
 );

 image_derivative_y(
  u,
  uy,
  xdim,
  ydim
 );

 image_derivative_x(
  v,
  vx,
  xdim,
  ydim
 );

 image_derivative_y(
  v,
  vy,
  xdim,
  ydim
 );

 /*
 Compute the derivative of psi for the smoothness term
 */

 psi_prime_smooth(
  ux,
  uy,
  vx,
  vy,
  psis,
  xdim,
  ydim,
  epsilon
 );

 /*
 Compute the terms needed to evaluate the divergence
 for the smoothness term
 */

 psi_prime_smooth_1234(
  psis,
  psi1,
  psi2,
  psi3,
  psi4,
  xdim,
  ydim
 );

 /*
 Compute the divergence of the gradient of the displacement
 */

 divergence_gradient_u(
  u,
  psi1,
  psi2,
  psi3,
  psi4,
  div_u,
  xdim,
  ydim
 );

 divergence_gradient_v(
  v,
  psi1,
  psi2,
  psi3,
  psi4,
  div_v,
  xdim,
  ydim
 );

 /*
 Compute the constant term in
 the divergence of the gradient of the displacement increment
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       div_d[i*xdim+j]=
        psi1[i*xdim+j]+
        psi2[i*xdim+j]+
        psi3[i*xdim+j]+
        psi4[i*xdim+j];
    }
 }

 /*
 Initialize the displacement increment
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       du[i*xdim+j]= 0;
       dv[i*xdim+j]= 0;
    }
 }

 inner_iter= 0;

 START2:

 /*
 Inner iteration
 */

 inner_iter++;

 /*
 Compute the robust version of the data term
 */

 psi_prime_data(
  I1,
  I2w,
  I2wx,
  I2wy,
  du,
  dv,
  psid,
  xdim,
  ydim,
  epsilon
 );

 /*
 Compute the robust version of the gradient term
 */

 psi_prime_gradient(
  I1,
  I1x,
  I1y,
  I2w,
  I2wx,
  I2wy,
  I2wxx,
  I2wxy,
  I2wyy,
  du,
  dv,
  psig,
  xdim,
  ydim,
  epsilon
 );

 /*
 Compute the terms in the du and dv equation
 */

 equation_du_dv(
  I1,
  I1x,
  I1y,
  I2w,
  I2wx,
  I2wy,
  I2wxx,
  I2wxy,
  I2wyy,
  psid,
  psig,
  div_u,
  div_v,
  div_d,
  Au,
  Av,
  Du,
  Dv,
  D,
  xdim,
  ydim,
  gamma,
  alpha
 );

 /*
 Compute du,dv using SOR iterations
 */

 sor_iteration(
  Au,
  Av,
  Du,
  Dv,
  D,
  du,
  dv,
  psi1,
  psi2,
  psi3,
  psi4,
  xdim,
  ydim,
  alpha,
  sor_w,
  max_sor_iter,
  max_sor_error
 );

 if ( inner_iter < max_inner_iter )
  goto START2;

 /*
 We are done with the inner iterations, meaning,
 we are done iterating on du,dv and we can update the displacement
 */

 /*
 Update the displacement
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       u[i*xdim+j]+= du[i*xdim+j];
       v[i*xdim+j]+= dv[i*xdim+j];
    }
 }

 /*
 Compute the min and max values for the displacement increment
 */

 min_du_val= 1.0e32;
 max_du_val=-1.0e32;
 avg_du_val=0.0;
 min_dv_val= 1.0e32;
 max_dv_val=-1.0e32;
 avg_dv_val=0.0;
 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       du_val= du[i*xdim+j];
       avg_du_val+= fabs(du_val);
       if ( du_val < min_du_val )
        min_du_val= du_val;
       if ( du_val > max_du_val )
        max_du_val= du_val;
       dv_val= dv[i*xdim+j];
       avg_dv_val=+ fabs(dv_val);
       if ( dv_val < min_dv_val )
        min_dv_val= dv_val;
       if ( dv_val > max_dv_val )
        max_dv_val= dv_val;
    }
 }
 avg_du_val/= (double)(xdim*ydim);
 avg_dv_val/= (double)(xdim*ydim);

 /* 
 fprintf(stdout,"min du = %8f max du = %8f avg du = %8f\n",min_du_val,max_du_val,avg_du_val);
 fprintf(stdout,"min dv = %8f max dv = %8f avg dv = %8f\n",min_dv_val,max_dv_val,avg_dv_val);
 */
 
 if ( outer_iter < max_outer_iter )
  goto START;

 /*
 We are done with the outer iterations, meaning,
 we are done warping the 2nd image with the displacement
 */

 /*
 Free all the allocated memory
 */

 free(I2w);
 free(I2wx);
 free(I2wy);
 free(I2wxx);
 free(I2wxy);
 free(I2wyy);
 free(ux);
 free(uy);
 free(vx);
 free(vy);
 free(psis);
 free(psi1);
 free(psi2);
 free(psi3);
 free(psi4);
 free(div_u);
 free(div_v);
 free(div_d);
 free(du);
 free(dv);
 free(psid);
 free(psig);
 free(Au);
 free(Av);
 free(Du);
 free(Dv);
 free(D);

}

#include "header.h"
#include "proto.h"

void compute_disparity_map(
 int *I1_int,
 int *I2_int,
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
 double max_sor_error,
 int **pdisp_u_arr,
 int *pmin_d_u,
 int *pmax_d_u,
 int **pdisp_v_arr,
 int *pmin_d_v,
 int *pmax_d_v
)

{

 int xdim;
 int ydim;
 double *I1;
 double *I2;
 int i;
 int j;
 int pixel;
 int r;
 int g;
 int b;
 double val;
 double *disp_u;
 double *disp_v;
 double min_disp_u;
 double max_disp_u;
 double min_disp_v;
 double max_disp_v;
 int *disp_u_arr;
 int *disp_v_arr;
 int val_int;
 int min_d_u;
 int max_d_u;
 int min_d_v;
 int max_d_v;

 xdim= w;
 ydim= h;

 /*
 Allocate memory to store image 1
 */

 I1= (double *)calloc(xdim*ydim,sizeof(double));

 /*
 Allocate memory to store image 2
 */

 I2= (double *)calloc(xdim*ydim,sizeof(double));

 /*
 Fill I1
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       pixel= i*xdim+j;
       r = I1_int[3*pixel+0];
       g = I1_int[3*pixel+1];
       b = I1_int[3*pixel+2];
       val= .2989*(double)r+.5870*(double)g+.1140*(double)b;
       I1[pixel]= val;
    }
 }

 /*
 Fill I2
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       pixel= i*xdim+j;
       r = I2_int[3*pixel+0];
       g = I2_int[3*pixel+1];
       b = I2_int[3*pixel+2];
       val= .2989*(double)r+.5870*(double)g+.1140*(double)b;
       I2[pixel]= val;
    }
 }

 disp_u= (double *)calloc(xdim*ydim,sizeof(double));
 disp_v= (double *)calloc(xdim*ydim,sizeof(double));

 brox_optic_flow_main(
  I1,
  I2,
  disp_u,
  disp_v,
  w,
  h,
  factor,
  nscales,
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
 Get the min and max u displacements
 */

 min_disp_u= 1.0e32;
 max_disp_u=-1.0e32;

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       val= disp_u[i*xdim+j];
       if ( val < min_disp_u )
        min_disp_u= val;
       if ( val > max_disp_u )
        max_disp_u= val;
    }
 }

 fprintf(stdout,"min_disp_u = %g max_disp_u = %g\n",min_disp_u,max_disp_u);

 /*
 Get the min and max v displacements
 */

 min_disp_v= 1.0e32;
 max_disp_v=-1.0e32;

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       val= disp_v[i*xdim+j];
      if ( val < min_disp_v )
        min_disp_v= val;
       if ( val > max_disp_v )
        max_disp_v= val;
    }
 }

 fprintf(stdout,"min_disp_v = %g max_disp_v = %g\n",min_disp_v,max_disp_v);

 /*
 Convert u displacement to disparity
 */

 disp_u_arr= (int *)calloc(xdim*ydim,sizeof(int));

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       pixel= i*xdim+j;
       val= disp_u[pixel];
       if ( val >= 0. )
        val_int= (int)(val+0.5);
       else
        val_int= (int)(val-0.5);
       disp_u_arr[pixel]= -val_int;
    }
 }

 /*
 Get the min and max disparities
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       pixel= i*xdim+j;
       val_int= disp_u_arr[pixel];
       if ( pixel == 0 ) {
          min_d_u= val_int;
          max_d_u= val_int;
       }
       else {
          if ( val_int < min_d_u )
           min_d_u= val_int;
          if ( val_int > max_d_u )
           max_d_u= val_int;
       }
    }
 }

 fprintf(stdout,"min_d_u = %d max_d_u = %d\n",min_d_u,max_d_u);

 /*
 Convert v displacement to disparity
 */

 disp_v_arr= (int *)calloc(xdim*ydim,sizeof(int));

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       pixel= i*xdim+j;
       val= disp_v[pixel];
       if ( val >= 0. )
        val_int= (int)(val+0.5);
       else
        val_int= (int)(val-0.5);
       disp_v_arr[pixel]= -val_int;
    }
 }

 /*
 Get the min and max disparities
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       pixel= i*xdim+j;
       val_int= disp_v_arr[pixel];
       if ( pixel == 0 ) {
          min_d_v= val_int;
          max_d_v= val_int;
       }
       else {
          if ( val_int < min_d_v )
           min_d_v= val_int;
          if ( val_int > max_d_v )
           max_d_v= val_int;
       }
    }
 }

 fprintf(stdout,"min_d_v = %d max_d_v = %d\n",min_d_v,max_d_v);

 /*
 Free memory
 */

 free(I1);
 free(I2);
 free(disp_u);
 free(disp_v);

 (*pdisp_u_arr)= disp_u_arr;
 (*pdisp_v_arr)= disp_v_arr;
 (*pmin_d_u)= min_d_u;
 (*pmax_d_u)= max_d_u;
 (*pmin_d_v)= min_d_v;
 (*pmax_d_v)= max_d_v;

}

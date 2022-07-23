#include "header.h"
#include "proto.h"

void gaussian_blur_image(
 double *I,
 int xdim,
 int ydim,
 double sigma,
 int precision
)

/*
The output image is put into the input image
*/

{

 int size;
 double *G;
 int i;
 double x;
 double pi=acos(-1);
 double *row;
 int j;
 double norm;
 int k;
 double val;
 double *col;

 /*
 We are gonna apply a 1d kernel
 (defined by sigma and precision)
 in both directions
 */

 size= (int)(precision*sigma)+1;

 G= (double *)calloc(size,sizeof(double));
 for ( i= 0 ; i< size ; i++ ) {
    x= (double)i;
    G[i]= 1/(sqrt(2*pi)*sigma)*exp(-x*x/(2*sigma*sigma));
 }

 /*
 We are gonna copy the image row by row,
 convolute each row, and
 put the results back in I
 */

 row= (double *)calloc(xdim,sizeof(double));
 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       norm= 0.0;
       k= 0;
       val= G[k]*I[i*xdim+j];
       norm+= G[k];
       for ( k= 1 ; k< size ; k++ ) {
          if ( !( (j+k)<xdim ) )
           continue;
          val+= G[k]*I[i*xdim+(j+k)];
          norm+= G[k];
       }
       for ( k= 1 ; k< size ; k++ ) {
          if ( !( (j-k)>=0 ) )
           continue;
          val+= G[k]*I[i*xdim+(j-k)];
          norm+= G[k];
       }
       row[j]= val/norm;
    }
    for ( j= 0 ; j< xdim ; j++ ) {
       I[i*xdim+j]= row[j];
    }
 }
 free(row);

 /*
 We are gonna copy the image col by col,
 convolute each col, and
 put the results back in I
 */

 col= (double *)calloc(ydim,sizeof(double));
 for ( j= 0 ; j< xdim ; j++ ) {
    for ( i= 0 ; i< ydim ; i++ ) {
       norm= 0.0;
       k= 0;
       val= G[k]*I[i*xdim+j];
       norm+= G[k];
       for ( k= 1 ; k< size ; k++ ) {
          if ( !( (i+k)<ydim ) )
           continue;
          val+= G[k]*I[(i+k)*xdim+j];
          norm+= G[k];
       }
       for ( k= 1 ; k< size ; k++ ) {
          if ( !( (i-k)>=0 ) )
           continue;
          val+= G[k]*I[(i-k)*xdim+j];
          norm+= G[k];
       }
       col[i]= val/norm;
    }
    for ( i= 0 ; i< ydim ; i++ ) {
       I[i*xdim+j]= col[i];
    }
 }
 free(col);

 free(G);

}

#include "header.h"
#include "proto.h"

void normalize_image_pair(
 double *inp_I1,
 double *inp_I2,
 double *out_I1,
 double *out_I2,
 int size
)

/*
out_I1,out_I2 must be allocated outside this function
*/

{

 double min1;
 double max1;
 double min2;
 double max2;
 int ind;
 double val;
 double min12;
 double max12;
 double diff;

 /*
 Compute the min and max values for both images
 */

 min1= 1.0e32;
 max1=-1.0e32;

 for ( ind= 0 ; ind< size ; ind++ ) {
    val= inp_I1[ind];
    if ( val < min1 )
     min1= val;
    if ( val > max1 )
     max1= val;
 }

 min2= 1.0e32;
 max2=-1.0e32;

 for ( ind= 0 ; ind< size ; ind++ ) {
    val= inp_I2[ind];
    if ( val < min2 )
     min2= val;
    if ( val > max2 )
     max2= val;
 }

 min12= min1;
 if ( min2 < min12 )
  min12= min2;

 max12= max1;
 if ( max2 > max12 )
  max12= max2;

 diff= max12-min12;

 /*
 Normalize both images between 0 and 255
 */

 for ( ind= 0 ; ind< size ; ind++ ) {
    val= inp_I1[ind];
    out_I1[ind]= 255.0*(val-min12)/diff;
 }

 for ( ind= 0 ; ind< size ; ind++ ) {
    val= inp_I2[ind];
    out_I2[ind]= 255.0*(val-min12)/diff;
 }

}

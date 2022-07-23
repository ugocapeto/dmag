#include "header.h"
#include "proto.h"

void normalize_image_dbl(
 double *inp_I,
 int *out_I,
 int size,
 double min_val,
 double max_val
)

/*
out_I must be allocated outside this function
*/

{

 int ind;
 double val;
 double ratio;

 /*
 Normalize the image values between 0 and 255
 */

 for ( ind= 0 ; ind< size ; ind++ ) {
    val= inp_I[ind];
    ratio= (val-min_val)/(max_val-min_val);
    out_I[ind]= (int)(255.*ratio);
 }

}

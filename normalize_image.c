#include "header.h"
#include "proto.h"

void normalize_image(
 int *inp_I,
 int *out_I,
 int size,
 int min_val,
 int max_val
)

/*
out_I must be allocated outside this function
*/

{

 int ind;
 int val;
 double ratio;

 /*
 Normalize the image values between 0 and 255
 */

 for ( ind= 0 ; ind< size ; ind++ ) {
    val= inp_I[ind];
    ratio= (double)(val-min_val)/(double)(max_val-min_val);
    out_I[ind]= (int)(255*ratio);
 }

}

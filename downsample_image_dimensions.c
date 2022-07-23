#include "header.h"
#include "proto.h"

void downsample_image_dimensions(
 int inp_xdim,
 int inp_ydim,
 int *pout_xdim,
 int *pout_ydim,
 double factor
)

/*
factor should be between 1 and 1
*/

{

 int out_xdim;
 int out_ydim;

 out_xdim= (int)( (double)inp_xdim*factor+.5 );
 out_ydim= (int)( (double)inp_ydim*factor+.5 );

 (*pout_xdim)= out_xdim;
 (*pout_ydim)= out_ydim;

}

#include "header.h"
#include "proto.h"

int main(
 int argc,
 char* argv[]
)

{

 FILE *fp;
 char filename[80];
 char filename1[80];
 char filename2[80];
 char filename_disp[80];
 int *I1_int;
 int *I2_int;
 int width;
 int height;
 int nscales;
 double factor;
 double epsilon;
 double alpha;
 int max_outer_iter;
 int max_inner_iter;
 double sor_w;
 int max_sor_iter;
 double max_sor_error;
 double gamma;
 int err_flag;
 int *disp_arr;
 int min_d;
 int max_d;
 int *disp_normalized_arr;

 /*
 Let's read in the input file
 */

 fp= fopen("dmag_input.txt","r");
 if ( fp == 0 ) {
    fprintf(stdout,"dmag_input.txt not found!\n");
    return 1;
 }

 /*
 Get filename for image 1
 */

 fscanf(fp,"%s",filename);

 fprintf(stdout,"image 1 = %s\n",filename);

 strcpy(filename1,filename);

 /*
 Get filename for image 2
 */

 fscanf(fp,"%s",filename);

 fprintf(stdout,"image 2 = %s\n",filename);

 strcpy(filename2,filename);

 /*
 Get the downsampling factor
 */

 fscanf(fp,"%lg",&factor);

 fprintf(stdout,"downsampling factor = %g\n",factor);

 /*
 Get the number of scales
 */

 fscanf(fp,"%d",&nscales);

 fprintf(stdout,"number of scales = %d\n",nscales);

 /*
 Get epsilon used in the robust terms
 */

 fscanf(fp,"%lg",&epsilon);

 fprintf(stdout,"epsilon (robust) = %g\n",epsilon);

 /*
 Get gamma (gradient term weight)
 */

 fscanf(fp,"%lg",&gamma);

 fprintf(stdout,"gamma (gradient) = %g\n",gamma);

 /*
 Get alpha (smoothness term weight)
 */

 fscanf(fp,"%lg",&alpha);

 fprintf(stdout,"alpha (smooth) = %g\n",alpha);

 /*
 Get max number of outer iterations
 */

 fscanf(fp,"%d",&max_outer_iter);

 fprintf(stdout,"max nbr of iterations (outer) = %d\n",max_outer_iter);

 /*
 Get max number of inner iterations
 */

 fscanf(fp,"%d",&max_inner_iter);

 fprintf(stdout,"max nbr of iterations (inner) = %d\n",max_inner_iter);

 /*
 Get sor relaxation factor
 */

 fscanf(fp,"%lg",&sor_w);

 fprintf(stdout,"relaxation factor = %g\n",sor_w);

 /*
 Get max number of sor iterations
 */

 fscanf(fp,"%d",&max_sor_iter);

 fprintf(stdout,"max nbr of iterations (sor) = %d\n",max_sor_iter);

 /*
 Get max sor error
 */

 fscanf(fp,"%lg",&max_sor_error);

 fprintf(stdout,"max sor error = %g\n",max_sor_error);

 /*
 Get filename for disparity map
 */

 fscanf(fp,"%s",filename_disp);

 fprintf(stdout,"disparity map = %s\n",filename_disp);

 /*
 We are done reading the input file and
 loading the 2 images
 */

 fclose(fp);

 /*
 Load image 1
 */

 err_flag= load_rgb_image(
  filename1,
  &I1_int,
  &width,
  &height
 );

 if ( err_flag == 1 )
  return 1;

 /*
 Load image 2
 */

 err_flag= load_rgb_image(
  filename2,
  &I2_int,
  &width,
  &height
 );

 if ( err_flag == 1 )
  return 1;

 /*
 Compute the number of scales if not given
 */

 if ( nscales == 0 ) {

    /*
    Let's compute the number of scales
    so that the image at scale (nscales-1)
    is at least 16 in width
    */

    nscales= 1+log((double)width/16.)/log(1./factor);

    fprintf(stdout,"number of scales = %d\n",nscales);
 }

 /*
 Ready to get the disparity maps
 */

 compute_disparity_map(
  I1_int,
  I2_int,
  width,
  height,
  factor,
  nscales,
  epsilon,
  gamma,
  alpha,
  max_inner_iter,
  max_outer_iter,
  sor_w,
  max_sor_iter,
  max_sor_error,
  &disp_arr,
  &min_d,
  &max_d
 );

 /*
 Let's dump the disparity map
 */

 disp_normalized_arr= (int *)calloc(width*height,sizeof(int));

 normalize_image(
  disp_arr,
  disp_normalized_arr,
  width*height,
  min_d,
  max_d
 );

 err_flag= write_image(
  filename_disp,
  disp_normalized_arr,
  width,
  height
 );

 if ( err_flag == 1 ) {
    return 1;
 }

 free(disp_normalized_arr);

 /*
 Free memory
 */

 free(I1_int);
 free(I2_int);
 free(disp_arr);

 return(0);

}

#include "header.h"
#include "proto.h"

void sor_iteration(
 double *Au,
 double *Av,
 double *Du,
 double *Dv,
 double *D,
 double *du,
 double *dv,
 double *psi1,
 double *psi2,
 double *psi3,
 double *psi4,
 int xdim,
 int ydim,
 double alpha,
 double sor_w,
 double max_sor_iter,
 double max_sor_error
)

/*
In the paper,
i is associated with x and j is associated with y
but in the code,
i is associated with y and j is associated with x
which means that
i must be swapped with j and j must be swapped with i
*/

{

 int sor_iter;
 double sum_val;
 double avg_val;
 int i;
 int j;
 double div_du;
 double div_dv;
 int ind;
 double valu;
 double valv;
 double val2;
 double val;

 sor_iter= 0;

 START:

 sor_iter++;
 sum_val= 0.0;

 for ( i= 0 ; i< ydim ; i++ ) { 
    for ( j= 0 ; j< xdim ; j++ ) {

       /*
       Compute div_du
       */

       div_du= 0.0;

       if ( (j+1) < xdim )
        div_du+= psi1[i*xdim+j]*du[i*xdim+(j+1)];
       else
        div_du+= psi1[i*xdim+j]*du[i*xdim+(j+0)];

       if ( (j-1) >= 0 )
        div_du+= psi2[i*xdim+j]*du[i*xdim+(j-1)];
       else
        div_du+= psi2[i*xdim+j]*du[i*xdim+(j-0)];

       if ( (i+1) < ydim )
        div_du+= psi3[i*xdim+j]*du[(i+1)*xdim+j];
       else
        div_du+= psi3[i*xdim+j]*du[(i+0)*xdim+j];

       if ( (i-1) >= 0 )
        div_du+= psi4[i*xdim+j]*du[(i-1)*xdim+j];
       else
        div_du+= psi4[i*xdim+j]*du[(i-0)*xdim+j];

       /*
       Compute div_dv
       */

       div_dv= 0.0;

       if ( (j+1) < xdim )
        div_dv+= psi1[i*xdim+j]*dv[i*xdim+(j+1)];
       else
        div_dv+= psi1[i*xdim+j]*dv[i*xdim+(j+0)];

       if ( (j-1) >= 0 )
        div_dv+= psi2[i*xdim+j]*dv[i*xdim+(j-1)];
       else
        div_dv+= psi2[i*xdim+j]*dv[i*xdim+(j-0)];

       if ( (i+1) < ydim ) 
        div_dv+= psi3[i*xdim+j]*dv[(i+1)*xdim+j];
       else 
        div_dv+= psi3[i*xdim+j]*dv[(i+0)*xdim+j];

       if ( (i-1) >= 0 )
        div_dv+= psi4[i*xdim+j]*dv[(i-1)*xdim+j];
       else 
        div_dv+= psi4[i*xdim+j]*dv[(i-0)*xdim+j];

       /*
       Update the displacement increment
       */

       ind= i*xdim+j;

       /*
       Keep track of the old values
       */

       valu= du[ind];
       valv= dv[ind];

       du[ind]= (1.0-sor_w)*du[ind]+
        sor_w*(Au[ind]-D[ind]*dv[ind]+alpha*div_du)/Du[ind];
       dv[ind]= (1.0-sor_w)*dv[ind]+
        sor_w*(Av[ind]-D[ind]*du[ind]+alpha*div_dv)/Dv[ind];

       /*
       Keep track of how much the increments change
       */

       val2= 0.0;
       val2+= (valu-du[ind])*(valu-du[ind]);
       val2+= (valv-dv[ind])*(valv-dv[ind]);
       val= sqrt(val2);

       sum_val+= val;
    }
 }

 avg_val= sum_val/(double)(xdim*ydim);

 /*
 Check for termination
 */

 if ( avg_val > max_sor_error && sor_iter < max_sor_iter )
  goto START;

}

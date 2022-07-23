#include "header.h"
#include "proto.h"

void equation_du_dv(
 double *I1,
 double *I1x,
 double *I1y,
 double *I2w,
 double *I2wx,
 double *I2wy,
 double *I2wxx,
 double *I2wxy,
 double *I2wyy,
 double *psid,
 double *psig,
 double *div_u,
 double *div_v,
 double *div_d,
 double *Au,
 double *Av,
 double *Du,
 double *Dv,
 double *D,
 int xdim,
 int ydim,
 double gamma,
 double alpha
)

{

 int i;
 int j;
 int ind;

 for ( i= 0 ; i< ydim ; i++ ) { 
    for ( j= 0 ; j< xdim ; j++ ) {

       ind= i*xdim+j;

       Au[ind]= 0.0;
       Au[ind]-= (1.0-gamma)*psid[ind]*(I2w[ind]-I1[ind])*I2wx[ind];
       Au[ind]-= gamma*psig[ind]*(I2wx[ind]-I1x[ind])*I2wxx[ind];
       Au[ind]-= gamma*psig[ind]*(I2wy[ind]-I1y[ind])*I2wxy[ind];
       Au[ind]+= alpha*div_u[ind];

       Av[ind]= 0.0;
       Av[ind]-= (1.0-gamma)*psid[ind]*(I2w[ind]-I1[ind])*I2wy[ind];
       Av[ind]-= gamma*psig[ind]*(I2wx[ind]-I1x[ind])*I2wxy[ind];
       Av[ind]-= gamma*psig[ind]*(I2wy[ind]-I1y[ind])*I2wyy[ind];
       Av[ind]+= alpha*div_v[ind];

       Du[ind]= 0.0;
       Du[ind]+= (1.0-gamma)*psid[ind]*I2wx[ind]*I2wx[ind];
       Du[ind]+= gamma*psig[ind]*I2wxx[ind]*I2wxx[ind];
       Du[ind]+= gamma*psig[ind]*I2wxy[ind]*I2wxy[ind];
       Du[ind]+= alpha*div_d[ind];

       Dv[ind]= 0.0;
       Dv[ind]+= (1.0-gamma)*psid[ind]*I2wy[ind]*I2wy[ind];
       Dv[ind]+= gamma*psig[ind]*I2wyy[ind]*I2wyy[ind];
       Dv[ind]+= gamma*psig[ind]*I2wxy[ind]*I2wxy[ind];
       Dv[ind]+= alpha*div_d[ind];

       D[ind]= 0.0;
       D[ind]+= (1.0-gamma)*psid[ind]*I2wx[ind]*I2wy[ind];
       D[ind]+= gamma*psig[ind]*I2wxx[ind]*I2wxy[ind];
       D[ind]+= gamma*psig[ind]*I2wyy[ind]*I2wxy[ind];
    }
 }

}

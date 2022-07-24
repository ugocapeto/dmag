void normalize_image_pair(
 double *inp_I1,
 double *inp_I2,
 double *out_I1,
 double *out_I2,
 int size
);

void gaussian_blur_image(
 double *I,
 int xdim,
 int ydim,
 double sigma,
 int precision
);

void downsample_image_dimensions(
 int inp_xdim,
 int inp_ydim,
 int *pout_xdim,
 int *pout_ydim,
 double factor
);

void downsample_image(
 double *inp_I,
 int inp_xdim,
 int inp_ydim,
 double *out_I,
 int out_xdim,
 int out_ydim,
 double factor
);

void warp_image(
 double *inp_I,
 double *u,
 double *v,
 double *out_I,
 int xdim,
 int ydim
);

void image_derivative_x(
 double *I,
 double *Ix,
 int xdim,
 int ydim
);

void image_derivative_y(
 double *I,
 double *Iy,
 int xdim,
 int ydim
);

void image_derivative_xx(
 double *I,
 double *Ixx,
 int xdim,
 int ydim
);

void image_derivative_xy(
 double *I,
 double *Ixy,
 int xdim,
 int ydim
);

void image_derivative_yy(
 double *I,
 double *Iyy,
 int xdim,
 int ydim
);

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
);

void upsample_image(
 double *inp_I,
 int inp_xdim,
 int inp_ydim,
 double *out_I,
 int out_xdim,
 int out_ydim
);

void psi_prime_smooth(
 double *ux,
 double *uy,
 double *vx,
 double *vy,
 double *psis,
 int xdim,
 int ydim,
 double epsilon
);

void psi_prime_smooth_1234(
 double *psis,
 double *psi1,
 double *psi2,
 double *psi3,
 double *psi4,
 int xdim,
 int ydim
);

void divergence_gradient_u(
 double *u,
 double *psi1,
 double *psi2,
 double *psi3,
 double *psi4,
 double *div_u,
 int xdim,
 int ydim
);

void divergence_gradient_v(
 double *v,
 double *psi1,
 double *psi2,
 double *psi3,
 double *psi4,
 double *div_v,
 int xdim,
 int ydim
);

void psi_prime_data(
 double *I1,
 double *I2w,
 double *I2wx,
 double *I2wy,
 double *du,
 double *dv,
 double *psid,
 int xdim,
 int ydim,
 double epsilon
);

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
);

void normalize_image(
 int *inp_I,
 int *out_I,
 int size,
 int min_val,
 int max_val
);

void normalize_image_dbl(
 double *inp_I,
 int *out_I,
 int size,
 double min_val,
 double max_val
);

void brox_optic_flow(
 double *I1,
 double *I2,
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
);

void brox_optic_flow_main(
 double *I1,
 double *I2,
 double *disp_u,
 double *disp_v,
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
 double max_sor_error
);

void psi_prime_gradient(
 double *I1,
 double *I1x,
 double *I1y,
 double *I2w,
 double *I2wx,
 double *I2wy,
 double *I2wxx,
 double *I2wxy,
 double *I2wyy,
 double *du,
 double *dv,
 double *psig,
 int xdim,
 int ydim,
 double epsilon
);

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
);

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
);

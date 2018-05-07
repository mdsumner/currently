#include <Rcpp.h>
using namespace Rcpp;

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Utils.h>


// TODO: imin2,imax2,ftrunc are replaced with Rf_ versions here
// is that we we should do?


// One advection time step using a 4th order Runge-Kutta scheme
//void runge_kutta(double*, double*, double*, double*, int*, int*, double*);

// http://web.cs.ucdavis.edu/~ma/ECS177/papers/particle_tracing.pdf
//
// Advance a single particle by one step using a 4th order runge kutta algorithm
//void runge_kutta_cpp(double *x, double *y, double *U, double *V, int *m, int *n, double *t_step){
// [[Rcpp::export]]
List rk(NumericVector x, NumericVector y, NumericVector U, NumericVector V,
        IntegerVector m, IntegerVector n, NumericVector t_step) {

  NumericVector xx(x.length());
  NumericVector yy(y.length());

  int lower_x_coord,lower_y_coord;
  double dx, dy;
  double u_velo1, u_velo2, u_velo3, u_velo4;
  double v_velo1, v_velo2, v_velo3, v_velo4;
  double u_velo, v_velo;
  double tmp_x, tmp_y;
  int k;

  for (int iparticle = 0; iparticle < x.length(); iparticle++) {
    //-----------------------------------------------
    // First coefficient
    // Get coordinates
    double coord_x = x[iparticle];
    double coord_y = y[iparticle];
    lower_x_coord = (int) Rf_ftrunc(coord_x);
    lower_y_coord = (int) Rf_ftrunc(coord_y);

    // Get distances and therefore weights
    dx = coord_x - lower_x_coord;
    dy = coord_y - lower_y_coord;

    // Take care of exceptions!
    lower_x_coord = Rf_imin2(Rf_imax2(lower_x_coord,1),(n[0] -1));
    lower_y_coord = Rf_imin2(Rf_imax2(lower_y_coord,1),(m[0] -1));

    // Get velocities
    k = (lower_x_coord -1)*(m[0]) + lower_y_coord -1;
    // -1 because array in numbered 0..n-1
    u_velo1 = (1-dx)*(1-dy)*U[k] + (1-dy)*dx*U[k+(m[0])] + dy*(1-dx)*U[k+1] + dx*dy*U[k+(m[0])+1];
    v_velo1 = (1-dx)*(1-dy)*V[k] + (1-dy)*dx*V[k+(m[0])] + dy*(1-dx)*V[k+1] + dx*dy*V[k+(m[0])+1];

    //-----------------------------------------------
    // Second coefficient
    // Get coordinates
    tmp_x = coord_x + (t_step[0] * u_velo1 / 2);
    tmp_y = coord_y + (t_step[0] * v_velo1 / 2);
    lower_x_coord = (int) Rf_ftrunc(tmp_x);
    lower_y_coord = (int) Rf_ftrunc(tmp_y);

    // Get distances and therefore weights
    dx = tmp_x - lower_x_coord;
    dy = tmp_y - lower_y_coord;

    // Take care of exceptions!
    lower_x_coord = Rf_imin2(Rf_imax2(lower_x_coord,1),(n[0] -1));
    lower_y_coord = Rf_imin2(Rf_imax2(lower_y_coord,1),(m[0] -1));

    // Get velocities
    k = (lower_x_coord -1)*(m[0]) + lower_y_coord -1;
    // -1 because array in numbered 0..n-1
    u_velo2 = (1-dx)*(1-dy)*U[k] + (1-dy)*dx*U[k+(m[0])] + dy*(1-dx)*U[k+1] + dx*dy*U[k+(m[0])+1];
    v_velo2 = (1-dx)*(1-dy)*V[k] + (1-dy)*dx*V[k+(m[0])] + dy*(1-dx)*V[k+1] + dx*dy*V[k+(m[0])+1];

    //-----------------------------------------------
    // Third coefficient
    // Get coordinates
    tmp_x = coord_x + (t_step[0] * u_velo2 / 2);
    tmp_y = coord_y + (t_step[0] * v_velo2 / 2);
    lower_x_coord = (int) Rf_ftrunc(tmp_x);
    lower_y_coord = (int) Rf_ftrunc(tmp_y);

    // Get distances and therefore weights
    dx = tmp_x - lower_x_coord;
    dy = tmp_y - lower_y_coord;

    // Take care of exceptions!
    lower_x_coord = Rf_imin2(Rf_imax2(lower_x_coord,1),(n[0] -1));
    lower_y_coord = Rf_imin2(Rf_imax2(lower_y_coord,1),(m[0] -1));

    // Get velocities
    k = (lower_x_coord -1)*(m[0]) + lower_y_coord -1;
    // -1 because array in numbered 0..n-1
    u_velo3 = (1-dx)*(1-dy)*U[k] + (1-dy)*dx*U[k+(m[0])] + dy*(1-dx)*U[k+1] + dx*dy*U[k+(m[0])+1];
    v_velo3 = (1-dx)*(1-dy)*V[k] + (1-dy)*dx*V[k+(m[0])]+ dy*(1-dx)*V[k+1] + dx*dy*V[k+(m[0])+1];

    //-----------------------------------------------
    // Fourth coefficient
    // Get coordinates
    tmp_x = coord_x + (t_step[0] * u_velo3);
    tmp_y = coord_y + (t_step[0] * v_velo3);
    lower_x_coord = (int) Rf_ftrunc(tmp_x);
    lower_y_coord = (int) Rf_ftrunc(tmp_y);

    // Get distances and therefore weights
    dx = tmp_x - lower_x_coord;
    dy = tmp_y - lower_y_coord;

    // Take care of exceptions!
    lower_x_coord = Rf_imin2(Rf_imax2(lower_x_coord,1),(n[0] -1));
    lower_y_coord = Rf_imin2(Rf_imax2(lower_y_coord,1),(m[0] -1));

    // Get velocities
    k = (lower_x_coord -1)*(m[0]) + lower_y_coord -1;
    // -1 because array in numbered 0..n-1
    u_velo4 = (1-dx)*(1-dy)*U[k] + (1-dy)*dx*U[k+(m[0])] + dy*(1-dx)*U[k+1] + dx*dy*U[k+(m[0])+1];
    v_velo4 = (1-dx)*(1-dy)*V[k] + (1-dy)*dx*V[k+(m[0])] + dy*(1-dx)*V[k+1] + dx*dy*V[k+(m[0])+1];


    //-----------------------------------------------
    // Combine coefficients
    u_velo = u_velo1/6 + u_velo2/3 + u_velo3/3 + u_velo4/6;
    v_velo = v_velo1/6 + v_velo2/3 + v_velo3/3 + v_velo4/6;

    // Move particles, add some random diffusion

    xx[iparticle] = coord_x + (u_velo * t_step[0]);
    yy[iparticle] = coord_y + (v_velo * t_step[0]);
  }


  return Rcpp::List::create(Rcpp::Named("x") = xx,
                            Rcpp::Named("y") = yy);
}

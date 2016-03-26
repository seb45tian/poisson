#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Convergence criteria */
#define eps 1e-06
#define kappa 100

/* POSITIONS OF THE POSITIVE AND NEGATIVE POLE */
double xpos = 0.3;
double ypos = 0.3;
double xneg = 0.7;
double yneg = 0.7;

/* Calculate kappa/pi here once (save some computation time) */
const double kappa_div_pi = kappa/M_PI;

/*=======================================================================*/
/* FUNCTION TO CALCULATE THE SQUARED MAGNITUDE OF TWO VECTORS: |r1-r2|^2 */
/*=======================================================================*/
double magnitude_squared(double x1, double y1, double x2, double y2) {

    double xret = x1-x2;
    double yret = y1-y2; 
    return((xret*xret)+(yret*yret));
}

/*===========================================================================*/
/* FUNCTION TO CALCULATE THE CORRECT CHARGE DENSITY RHO FOR A GIVEN VECTOR R */
/*===========================================================================*/
double rho_val(double x, double y) {
    double ret = kappa_div_pi * ( exp(-kappa*magnitude_squared(x,y,xpos,ypos)) 
                                    - exp(-kappa*magnitude_squared(x,y,xneg,yneg)) );
    return(ret);
}

/*=============================================================*/
/* FUNCTION TO CALCULATE THE MAGNITUDE OF TWO VECTORS: |r1-r2| */
/*=============================================================*/
double magnitude(double x1, double y1, double x2, double y2) {
    double xret = x1-x2;
    double yret = y1-y2;
    double sol  = sqrt((xret*xret)+(yret*yret));
    return(sol);
}

/*==============================================*/
/* FUNCTION TO CALCULATE MAX NORM: max(|r1-r2|) */
/*==============================================*/
double get_max(double **a, double **b, int Nx, int Ny) {
    double max = 0.0;
    for (int i=0;i<Nx;i++) {
        for (int j=0;j<Ny;j++)
            if (fabs(a[i][j]-b[i][j]) > max) {
                max =  fabs(a[i][j]-b[i][j]);
            }
    }
    return(max);
}


/*=======================================================*/
/* ITERATIVE POISSON FUNCTION TO CALCULATE THE POTENTIAL */
/*=======================================================*/
void poisson(double **phi, double **phi_new, double **phi_swap, double h, int N) {
    int i,j,iter;
    double max_norm;
    int iter_max = 1e05;

    for(iter=0; iter<iter_max; iter++) {
        for (i=0;i<N;i++) {
            for (j=0;j<N;j++) {
                /* the normal case when at no boundary */
                if (i!=0 && i!=N-1 && j!=0 && j!=N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - h*h*rho_val(h*i, h*j));
                }
                /* when at left-top corner */
                else if (i==0 && j==0) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1+N][j] + phi[i][j+1] + phi[i][j-1+N] - h*h*rho_val(h*i, h*j));
                }
                /* when at left but not at top or bottom boundary */
                else if (i==0 && j!=0 && j!=N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1+N][j] + phi[i][j+1] + phi[i][j-1] - h*h*rho_val(h*i, h*j)); 
                }
                /* when at top but not at left or right boundary */
                else if (i!=0 && j==0 && i!=N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1+N] - h*h*rho_val(h*i, h*j));
                }
                /* when at right-bottom corner */
                else if (i==N-1 && j==N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1-N][j] + phi[i-1][j] + phi[i][j+1-N] + phi[i][j-1] - h*h*rho_val(h*i, h*j));
                }
                /* when at right but not at top or bottom boundary */
                else if (i==N-1 && j!=N-1 && j!=0) {
                    phi_new[i][j] = 0.25*(phi[i+1-N][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - h*h*rho_val(h*i, h*j));
                }
                /* when at bottom but not at left or rigth boundary */
                else if (i!=N-1 && j==N-1 && i!=0) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1-N] + phi[i][j-1] - h*h*rho_val(h*i, h*j));
                }
                /* when at right-top corner */
                else if (i==N-1 && j==0) {
                    phi_new[i][j] = 0.25*(phi[i+1-N][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1+N] - h*h*rho_val(h*i, h*j));
                }
                /* when at left-bottom corner */
                else if (i==N-1 && j==0) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1+N][j] + phi[i][j+1-N] + phi[i][j-1] - h*h*rho_val(h*i, h*j));
                }

            }
        }
        /* Break loop if convergence reached */
        max_norm = get_max(phi,phi_new,N,N);
        if (max_norm <= eps) {
            printf("# Iterations: %d\n", iter);
            break;
        }

        /* Don't copy data - just swap pointers around */
        phi_swap = phi;
        phi = phi_new;
        phi_new = phi_swap;
    } /*end iter*/
}





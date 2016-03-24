#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "alloc.h"

/* Convergence criteria */
#define eps 1e-06
#define kappa 100

double xpos = 0.3;
double ypos = 0.3;
double xneg = 0.7;
double yneg = 0.7;


double magnitude_squared(double x1, double y1, double x2, double y2);
double magnitude(double x1, double y1, double x2, double y2);
double rho_val(double x, double y);
double get_max(double **a, int Nx, int Ny);


int main(int argc, char const *argv[])
{
    /* VARIBALES */
    int N = 5;             // number of cells
    double h = 1.0/N;       // grid spacing
    N++;                    // increase cells +1 to really get to the edges
    int xi = 0;             // initial value of xi
    int yj = 0;             // initial value of yi
    double xlength = 1.0;   // side length in x direction
    double ylength = 1.0;   // side length in y direction
    /* Create arrays */
    double **phi, **phi_new, **phi_swap,**rho;
    int i,j;

    /* ALLOCATE MEMORY */
    phi     = alloc_doublematrix(N, N);
    phi_new = alloc_doublematrix(N, N);
    rho     = alloc_doublematrix(N, N);

    /* Set all values to zero */
    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            phi[i][j] = 0.0;
            phi_new[i][j] = 0.0;
        }
    }


    for(int iter=0; iter<100; iter++) {
        for (i=0;i<N;i++) {
            for (j=0;j<N;j++) {
                if (i!=0 && i!=N-1 && j!=0 && j!=N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - h*h*rho_val(h*i, h*j));
                }
                else if (i==0 && j==0) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1+N][j] + phi[i][j+1] + phi[i][j-1+N] - h*h*rho_val(h*i, h*j));
                }
                else if (i==0 && j!=0 && j!=N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1+N][j] + phi[i][j+1] + phi[i][j-1] - h*h*rho_val(h*i, h*j)); 
                }
                else if (i!=0 && j==0 && i!=N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1+N] - h*h*rho_val(h*i, h*j));
                }
                else if (i==N-1 && j==N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1-N][j] + phi[i-1][j] + phi[i][j+1-N] + phi[i][j-1] - h*h*rho_val(h*i, h*j));
                }
                else if (i==N-1 && j!=N-1 && j!=0) {
                    phi_new[i][j] = 0.25*(phi[i+1-N][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - h*h*rho_val(h*i, h*j));
                }
                else if (i!=N-1 && j==N-1 && i!=0) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1-N] + phi[i][j-1] - h*h*rho_val(h*i, h*j));
                }

            }
        }
        /* Don't copy data - just swap pointers around */
        phi_swap = phi;
        phi = phi_new;
        phi_new = phi_swap;
    } /*end iter*/


    // for(i=0;i<N;i++) {
    //     for(j=0;j<N;j++) {
    //         printf("%f\n", phi[i][j]);
    //     }
    // }

    /* FREE MEMORY */
    free_matrix(phi);
    free_matrix(phi_new);
    free_matrix(rho);

    return 0;
}

double magnitude_squared(double x1, double y1, double x2, double y2) {

    double xret = x1-x2;
    double yret = y1-y2; 
    return((xret*xret)+(yret*yret));
}

double rho_val(double x, double y) {
    double rho_test = (kappa/M_PI) * ( exp(-kappa*magnitude_squared(x,y,xpos,ypos)) - exp(-kappa*magnitude_squared(x,y,xneg,yneg)) );
}

double magnitude(double x1, double y1, double x2, double y2) {
    double xret = x1-x2;
    double yret = y1-y2;
    double sol  = sqrt((xret*xret)+(yret*yret));
    return(sol);
}

double get_max(double **a, int Nx, int Ny) {
    double max = 0.0;
    for (int i=0;i<Nx;i++) {
        for (int j=0;j<Ny;j++)
            if (fabs(a[i][j]) > max) {
                max =  fabs(a[i][j]);
            }
    }
    return(max);
}







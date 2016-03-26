#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "alloc.h"
#include "simulation.h"


int main(int argc, char const *argv[])
{
    /* VARIBALES */
    int N = 51;                 // number of cells
    double length = 1.0;        // side length of the domain
    double h = length/(N-1);    // grid spacing
    /* Create data arrays */
    double **phi, **phi_new, **rho, **Ex, **Ey;

    /* ALLOCATE MEMORY (it's automatically allocated to 0 using calloc)*/
    phi     = alloc_doublematrix(N, N);
    phi_new = alloc_doublematrix(N, N);
    rho     = alloc_doublematrix(N, N);
    Ex      = alloc_doublematrix(N, N);
    Ey      = alloc_doublematrix(N, N);

    /* CALCULATE THE CHARGE DENSITY RHO OF THE DOMAIN */
    calc_density(rho,h,N);

    /* RUN THE POISSON CALCULATION */
    poisson(phi, phi_new, rho, h, N);

    /* CALCULATE THE FIELD E */
    calc_field(Ex, Ey, phi, h, N);
    
    print_matrix(phi, N, N, "Phi");
    // print_matrix(Ex, N, N, "Ex");
    // print_matrix(Ey, N, N, "Ey");

    /* FREE MEMORY */
    free_matrix(phi);
    free_matrix(phi_new);
    free_matrix(rho);
    free_matrix(Ex);
    free_matrix(Ey);

    return 0;
}









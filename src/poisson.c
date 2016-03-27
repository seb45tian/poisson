#include <stdio.h>
#include "alloc.h"
#include "simulation.h"


int main(int argc, char const *argv[])
{
    /* VARIBALES */
    int N = 51;                 // number of cells
    int iter_max = 1e5;         // number of maximum iterations
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
    poisson(phi, phi_new, rho, h, N, iter_max);

    /* CALCULATE THE FIELD E */
    calc_field(Ex, Ey, phi, h, N);
    
    /* WRITE DATA TO FILES */
    write_matrix(phi, N, "Phi", "potential.dat", "w");
    write_matrix(Ex, N, "Ex", "field.dat", "w");
    write_matrix(Ey, N, "Ey", "field.dat", "a");


    /* FREE MEMORY */
    free_matrix(phi);
    free_matrix(phi_new);
    free_matrix(rho);
    free_matrix(Ex);
    free_matrix(Ey);

    return 0;
}









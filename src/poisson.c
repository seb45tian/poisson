#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "alloc.h"
#include "simulation.h"


int main(int argc, char const *argv[])
{
    /* VARIBALES */
    int N = 50;             // number of cells
    double length = 1.0;    // side length of the domain
    double h = length/N;    // grid spacing
    N++;                    // increase cells +1 to really get to the edges
    /* Create arrays */
    double **phi, **phi_new, **phi_swap, **rho;
    int i,j;

    /* ALLOCATE MEMORY (it's automatically allocated to 0 using calloc)*/
    phi     = alloc_doublematrix(N, N);
    phi_new = alloc_doublematrix(N, N);
    rho     = alloc_doublematrix(N, N);

    /* RUN THE POISSON CALCULATION */
    poisson(phi, phi_new, phi_swap, h, N);


    print_matrix(phi, N, N, "Phi");

    /* FREE MEMORY */
    free_matrix(phi);
    free_matrix(phi_new);
    free_matrix(rho);

    return 0;
}









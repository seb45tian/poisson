#include <stdio.h>
#include <mpi.h>
#include "alloc.h"
#include "simulation.h"

int proc;                   /* rank of the current process */
int nprocs;                 /* total number of processes */
int *itop, *ibottom;        /* pointers for the boundariy values */

int main(int argc, char **argv)
{
    /* VARIBALES */
    int N = 51;                 // number of cells
    int iter_max = 1e4;         // number of maximum iterations
    double length = 1.0;        // side length of the domain
    double h = length/(N-1);    // grid spacing

    /* Create data arrays */
    double **phi, **phi_new, **rho, **Ex, **Ey;

    /*===================================================*/
    /* ADDED VARIABLES FOR MPI */
    int nrows;          /* The number of coulmns */
    int rem;            /* remainder of integer division imax/nprocs */
    int tpoint = 0;     /* stores the left column index for each proc*/
    int bpoint = 0;     /* stores the right column index for each proc*/
    double t1, t2;

    /* Initialization of pointers itop and ibottom */
    itop = &tpoint;
    ibottom = &bpoint; 

    /*===================================================*/
    /* INIT MPI*/
    /*===================================================*/
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &proc );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    t1 = MPI_Wtime(); // start the 
    /*===================================================*/



    /* ALLOCATE MEMORY (it's automatically allocated to 0 using calloc)*/
    phi     = alloc_doublematrix(N, N);
    phi_new = alloc_doublematrix(N, N);
    rho     = alloc_doublematrix(N, N);
    Ex      = alloc_doublematrix(N, N);
    Ey      = alloc_doublematrix(N, N);



    /*===================================================*/
    /* 1D Decomposition into rows */
    nrows = (int) N/nprocs;  /* number of rows for each proc */
    rem = N%nprocs;          /* get the remainder of the division */

    /* Every proc <= (rem-1) will get an extra column */
    if( proc < rem ) {
        nrows++;
    }
 
    /* Find the left starting index i (colum) for each processor */
    /* and save it in the tpoint variable                        */ 
    if( proc < rem ){
        tpoint = (proc*nrows);// + 1;
    } else {
        tpoint = (proc*nrows + rem);// + 1;  
    }

    /* Find the right end index i (colum) for each processor */
    /* and save it in the bpoint variable                    */  
    bpoint = (tpoint + nrows);// - 1;
    /*===================================================*/



    /* CALCULATE THE CHARGE DENSITY RHO OF THE DOMAIN */
    calc_density(rho,h,N);

    /* RUN THE POISSON CALCULATION */
    poisson(phi, phi_new, rho, h, N, iter_max);

    /* CALCULATE THE FIELD E */
    if (proc==0) {
        calc_field(Ex, Ey, phi, h, N);
    }
    
    /* WRITE DATA TO FILES */
    if (proc==0) {
        write_matrix(phi, N, "Phi", "potential.dat", "w");
        write_matrix(Ex, N, "Ex", "field.dat", "w");
        write_matrix(Ey, N, "Ey", "field.dat", "a");
    }


    /* FREE MEMORY */
    free_matrix(phi);
    free_matrix(phi_new);
    free_matrix(rho);
    free_matrix(Ex);
    free_matrix(Ey);


    /*===================================================*/
    /* CALCULATE THE COMPUTATION TIME */
    if(proc == 0){
        t2 = MPI_Wtime();
        printf("Used processors: %d, Time: %4.4f seconds\n", nprocs, (t2-t1) );
    }
    /* SHUT DOWN MPI*/
    MPI_Finalize();
    /*===================================================*/

    return 0;
}









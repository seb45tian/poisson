#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "alloc.h"

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

extern int *itop, *ibottom; /* pointers to the left and right column index */
extern int nprocs, proc;    /* total number of processors and rank id */
extern int force_iterations;

/*=======================================================================*/
/* FUNCTION TO CALCULATE THE SQUARED MAGNITUDE OF TWO VECTORS: |r1-r2|^2 */
/*=======================================================================*/
double magnitude_squared(double x1, double y1, double x2, double y2) {
    double xret = x1-x2;
    double yret = y1-y2; 
    return((xret*xret)+(yret*yret));
}


/*=============================================================*/
/* FUNCTION TO CALCULATE THE CHARGE DENSITY FOR A GIVEN DOMAIN */
/*=============================================================*/
void calc_density(double **rho, double h, int N) {
    int i,j;
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++) {
            rho[i][j] = kappa_div_pi * ( exp(-kappa*magnitude_squared(h*i,h*j,xpos,ypos)) 
                                    - exp(-kappa*magnitude_squared(h*i,h*j,xneg,yneg)) );
        }
    }
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
double get_max(double **a, double **b, int imin, int imax, int jmin, int jmax) {
    int i,j;
    double max = 0.0;
    for (i=imin; i<imax; i++) {
        for (j=jmin; j<jmax; j++) {
            if (fabs(a[i][j]-b[i][j]) > max) {
                max =  fabs(a[i][j]-b[i][j]);
            }
        }
    }
    return(max);
}


/*=======================================================*/
/* ITERATIVE POISSON FUNCTION TO CALCULATE THE POTENTIAL */
/*=======================================================*/
void poisson(double **phi, double **phi_new, double **rho, double h, int N, int iter_max) {
    int i,j,iter,ip;
    double max_norm;
    double max_norm_total;
    double hsq = h*h;
    double **phi_swap;


    /*===================================================*/
    MPI_Status status;      /* MPI status */
    int nrows;              /* number of rows */
    int rem;                /* remainder of integer division */

    int nrows_array[nprocs];        /* array holding the nrows of each proc */
    int tpoint_array[nprocs];       /* array holding the tpoint of each proc */
    
    /* GET THE CORRECT NUMBER OF ROWS AND REMAINDER */
    nrows = N/nprocs;
    rem = N%nprocs;
    
    /* Find the top starting row index of each proc and the number   */
    /* of rows and store them in arrays tpoint_array and nrows_array */
    /* Every proc needs to have this array for later use!            */
    for(ip =0; ip < nprocs; ip++){
        if(ip < rem){
            nrows_array[ip]  = nrows+1;
            tpoint_array[ip] = ip*nrows_array[ip];
        } else {
            nrows_array[ip]  = nrows;
            tpoint_array[ip] = ip*nrows_array[ip] + rem;
        }
    }
    /*===================================================*/






    for(iter=0; iter<iter_max; iter++) {
        for (i=*itop;i<*ibottom;i++) {
            for (j=0;j<N;j++) {
                /* the normal case when at no boundary */
                if (i!=0 && i!=N-1 && j!=0 && j!=N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - hsq*rho[i][j]);
                }
                /* when at left-top corner */
                else if (i==0 && j==0) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1+N][j] + phi[i][j+1] + phi[i][j-1+N] - hsq*rho[i][j]);
                }
                /* when at left but not at top or bottom boundary */
                else if (i==0 && j!=0 && j!=N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1+N][j] + phi[i][j+1] + phi[i][j-1] - hsq*rho[i][j]); 
                }
                /* when at top but not at left or right boundary */
                else if (i!=0 && j==0 && i!=N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1+N] - hsq*rho[i][j]);
                }
                /* when at right-bottom corner */
                else if (i==N-1 && j==N-1) {
                    phi_new[i][j] = 0.25*(phi[i+1-N][j] + phi[i-1][j] + phi[i][j+1-N] + phi[i][j-1] - hsq*rho[i][j]);
                }
                /* when at right but not at top or bottom boundary */
                else if (i==N-1 && j!=N-1 && j!=0) {
                    phi_new[i][j] = 0.25*(phi[i+1-N][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] - hsq*rho[i][j]);
                }
                /* when at bottom but not at left or rigth boundary */
                else if (i!=N-1 && j==N-1 && i!=0) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1][j] + phi[i][j+1-N] + phi[i][j-1] - hsq*rho[i][j]);
                }
                /* when at right-top corner */
                else if (i==N-1 && j==0) {
                    phi_new[i][j] = 0.25*(phi[i+1-N][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1+N] - hsq*rho[i][j]);
                }
                /* when at left-bottom corner */
                else if (j==N-1 && i==0) {
                    phi_new[i][j] = 0.25*(phi[i+1][j] + phi[i-1+N][j] + phi[i][j+1-N] + phi[i][j-1] - hsq*rho[i][j]);
                }

            } /*end of j*/
        } /*end of i*/


        /*======================================================*/
        /* WE HAVE TO SEND AND RECEIVE THE CORRECT BOUNDARIES   */
        /* TO THE PROC ABOVE AND BENEATH FROM US (proc+1 or -1) */

        if (nprocs > 1) {
            /* If inner processor: send and receive top AND bottom boundary */
            if( proc > 0 && proc < nprocs-1){ 
                MPI_Sendrecv(phi_new[*ibottom-1], N, MPI_DOUBLE, proc+1, 666, phi_new[*ibottom], N,
                             MPI_DOUBLE, proc+1, 999, MPI_COMM_WORLD, &status);
                MPI_Sendrecv(phi_new[*itop], N, MPI_DOUBLE, proc-1, 999, phi_new[*itop-1], N,
                             MPI_DOUBLE, proc-1, 666, MPI_COMM_WORLD, &status);
            /* If top or bottom boundary processor: only send and receive bottom OR top boundary */
            } else {
                if (proc == 0) {
                    MPI_Sendrecv(phi_new[*ibottom-1], N, MPI_DOUBLE, proc+1, 666, phi_new[*ibottom], N,
                                 MPI_DOUBLE, proc+1, 999, MPI_COMM_WORLD, &status);
                    /* for the periodic boundary conditions */
                    MPI_Sendrecv(phi_new[*itop], N, MPI_DOUBLE, nprocs-1, 333, phi_new[N-1], N,
                                 MPI_DOUBLE, nprocs-1, 444, MPI_COMM_WORLD, &status);        
                } else {
                    MPI_Sendrecv(phi_new[*itop], N, MPI_DOUBLE, proc-1, 999, phi_new[*itop-1], N,
                                 MPI_DOUBLE, proc-1, 666,MPI_COMM_WORLD, &status);
                    /* for the periodic boundary conditions */ 
                    MPI_Sendrecv(phi_new[*ibottom-1], N, MPI_DOUBLE, 0, 444, phi_new[0], N,
                                 MPI_DOUBLE, 0, 333, MPI_COMM_WORLD, &status);     
                }
            } 
        }
        /*===================================================*/
        /* CHECKING FOR CONVERGENCE */
        /*(Only if force_iterations is disabled)*/
        if (force_iterations == 0) {
            /* Calculate the local max-norm value */
            max_norm = get_max(phi_new, phi, *itop, *ibottom, 0, N);

            /* Get the maximum value of all max-norms */
            MPI_Allreduce(&max_norm, &max_norm_total, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

            /* Break loop if convergence reached */
            if (max_norm_total <= eps) {
                if (proc == 0) printf("# Reached convergence after %d iteration.\n", iter);
                break;
            }
        }
        /*===================================================*/

        /* Don't copy data - just swap pointers around */
        phi_swap = phi;
        phi = phi_new;
        phi_new = phi_swap;

    } /*end iter*/

    /*===================================================*/
    /* NOW WE COLLECT EVERYTHING ON PROC 0 */
    if (nprocs > 1) {
        if ( proc != 0){
            /* If not proc 0 --> send data */
            MPI_Send( phi[*itop], nrows_array[proc]*N, MPI_DOUBLE, 0, 777, MPI_COMM_WORLD);
        } else if (nprocs>1){
            /* If proc 0 --> loop over all other procs and receive data */
            for (ip = 1; ip < nprocs; ip++){
                MPI_Recv(phi[tpoint_array[ip]], nrows_array[ip]*N, MPI_DOUBLE, ip, 777, MPI_COMM_WORLD, &status);
            }
        }
         
        /* BROADCAST PHI ARRAY TO ALL PROCESSES */
        MPI_Bcast( phi[0],N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        /* WAIT UNTIL EVERYONE FINISHED */
        MPI_Barrier(MPI_COMM_WORLD);
    }
    /*===================================================*/

}



/*=========================================================================*/
/* FUNCTION TO CALCULATE THE FIELD COMPONENTS Ex and Ey FOR A GIVEN DOMAIN */
/*=========================================================================*/
void calc_field(double **Ex, double **Ey, double **phi, double h, int N) {
    int i,j;
    const double div2h = 1.0/(2*h);
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++) {
            // normal case, inner points
            if (i!=0 && i!=N-1 && j!=0 && j!=N-1) {
                Ex[i][j] = (phi[i+1][j]-phi[i-1][j])*div2h;
                Ey[i][j] = (phi[i][j+1]-phi[i][j-1])*div2h;
            }
            /* when at top-left corner */
            else if (i==0 && j==0) {
                Ex[i][j] = (phi[i+1][j]-phi[i-1+N][j])*div2h;
                Ey[i][j] = (phi[i][j+1]-phi[i][j-1+N])*div2h;
            }
            /* when at top-right corner */
            else if (i==N-1 && j==0) {
                Ex[i][j] = (phi[i+1-N][j]-phi[i-1][j])*div2h;
                Ey[i][j] = (phi[i][j+1]-phi[i][j-1+N])*div2h;
            }
            /* when at bottom-left corner */
            else if (i==0 && j==N-1) {
                Ex[i][j] = (phi[i+1][j]-phi[i-1+N][j])*div2h;
                Ey[i][j] = (phi[i][j+1-N]-phi[i][j-1])*div2h;
            }
            /* when at bottom-right corner */
            else if (i==N-1 && j==N-1) {
                Ex[i][j] = (phi[i+1-N][j]-phi[i-1][j])*div2h;
                Ey[i][j] = (phi[i][j+1-N]-phi[i][j-1])*div2h;
            }
            /* when at left boundary but not top/bottom corner */
            else if (i==0 && j!=0 && j!=N-1) {
                Ex[i][j] = (phi[i+1][j]-phi[i-1+N][j])*div2h;
                Ey[i][j] = (phi[i][j+1]-phi[i][j-1])*div2h;
            }
            /* when at right boundary but not top/bottom corner */
            else if (i==N-1 && j!=0 && j!=N-1) {
                Ex[i][j] = (phi[i+1-N][j]-phi[i-1][j])*div2h;
                Ey[i][j] = (phi[i][j+1]-phi[i][j-1])*div2h;
            }
            /* when at top boundary but not left/right corner */
            else if (j==0 && i!=0 && i!=N-1) {
                Ex[i][j] = (phi[i+1][j]-phi[i-1][j])*div2h;
                Ey[i][j] = (phi[i][j+1]-phi[i][j-1+N])*div2h;
            }
            /* when at bottom boundary but not left/right corner */
            else if (j==N-1 && i!=0 && i!=N-1) {
                Ex[i][j] = (phi[i+1][j]-phi[i-1][j])*div2h;
                Ey[i][j] = (phi[i][j+1-N]-phi[i][j-1])*div2h;
            }
        }
    }
}





#include <stdlib.h>
#include <stdio.h>

/* Allocate memory for a rows*cols array of floats.
 * The elements within a column are contiguous in memory, and columns
 * themselves are also contiguous in memory.
 */
double **alloc_doublematrix(int cols, int rows)
{
    int i;
    double **m;
    if ((m = (double**) malloc(cols*sizeof(double*))) == NULL) {
        return NULL;
    }
    double *els = (double *) calloc(rows*cols, sizeof(double));
    if (els == NULL) {
        return NULL;
    } 
    for (i = 0; i < cols; i++) {
        m[i] = &els[rows * i];
    }
    return m;
} 


/* Free the memory of a matrix allocated with alloc_doublematrix*/
void free_matrix(void *m)
{
    void **els = (void **) m;
    free(els[0]);   /* Deallocate the block of array elements */
    free(m);        /* Deallocate the block of column pointers */
}


// Function to print out an array
void print_matrix(double **data, int imax, int jmax, char *str) {    
    printf("-- %s --\n", str);
    for (int i=0; i<imax; i++) {
        for (int j=0; j<jmax; j++) {
            printf("%f ", data[i][j]);
        }
        printf("\n");
    }
}
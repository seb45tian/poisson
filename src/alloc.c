#include <stdlib.h>
#include <stdio.h>

/*========================================================*/
/* Allocate memory for a rows*cols array of floats.       */
/* The elements within a column are contiguous in memory, */ 
/* and columns themselves are also contiguous in memory.  */
/*========================================================*/
double **alloc_doublematrix(int cols, int rows) {
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

/*===============================================================*/
/* Free the memory of a matrix allocated with alloc_doublematrix */
/*===============================================================*/
void free_matrix(void *m) {
    void **els = (void **) m;
    free(els[0]);   /* Deallocate the block of array elements */
    free(m);        /* Deallocate the block of column pointers */
}

/*=======================================*/
/* Function to print 2D matrix to stdout */
/*=======================================*/
void print_matrix(double **data, int imax, int jmax, char *str) {  
    int i,j;  
    printf("#-- %s --#\n", str);
    for (i=0; i<imax; i++) {
        for (j=0; j<jmax; j++) {
            printf("%f ", data[i][j]);
        }
        printf("\n");
    }
}

/*======================================*/
/* Function to save 2D matrix to a file */
/*======================================*/
void write_matrix(double **data, int N, char *head, char *filename, char *mode) {
    int i,j;
    FILE *pFile;
    pFile = fopen(filename, mode);
    fprintf(pFile, "#-- %s --#\n", head);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            fprintf(pFile, "%f ", data[i][j]);
        }
        fprintf(pFile,"\n");
    }
    fclose (pFile);
}



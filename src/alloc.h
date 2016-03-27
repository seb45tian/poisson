double **alloc_doublematrix(int cols, int rows);
void free_matrix(void *m);
void print_matrix(double **data, int imax, int jmax, char *str);
void write_matrix(double **data, int N, char *head, char *filename, char *mode);
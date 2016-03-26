double magnitude_squared(double x1, double y1, double x2, double y2);
double magnitude(double x1, double y1, double x2, double y2);
double rho_val(double x, double y);
double get_max(double **a, double **b, int Nx, int Ny);
void poisson(double **phi, double **phi_new, double **phi_swap, double h, int N);
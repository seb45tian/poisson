double magnitude_squared(double x1, double y1, double x2, double y2);
double magnitude(double x1, double y1, double x2, double y2);
double get_max(double **a, double **b, int Nx, int Ny);
void calc_density(double **rho, double h, int N);
void poisson(double **phi, double **phi_new, double **rho, double h, int N);
void calc_field(double **Ex, double **Ey, double **phi, double h, int N);
void canonical2emittance_(double cancord[6], double emittance[3]);
void six2canonical_(double * coord, double *ref_momentum, double *mass, double *canonical);
void canonical2six_(double *canonical, double *ref_momentum, double *mass, double *coord);
double rationalApproximation(double t);
double normalcdfinv_(double p);
void createLinearSpaced(int length, double start, double stop, double *eqspaced );
double momentum2energy(double momentum, double mass);
void mtrx_vector_mult_(int *mp, int *np,  double mtrx_a[6][6], double mtrx_b[6], double result[6]);
void mtrx_vector_mult_pointer(int mp, int np,  double **mtrx_a, double mtrx_b[6], double result[6]);
void transpose(double num[6][6],double fac[6][6],double  r);
void cofactor(double num[6][6],double f);
void printmatrix(int m, int n, double **matrix );
void printvector(const char* name, int dim, double* vector);
void bisection(double *x, double a, double b, int *itr);
void hello();
double randn(double mu, double sigma);
double rand_uni();
void calcualteinverse();
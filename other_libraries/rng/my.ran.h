/* prototypes for my random number stuff */

extern double lchoose(double, double);
extern float lchoosef(float, float);

extern double runif();

extern double rgamma(double, double);
extern double dgamma(double, double, double, int);

extern void rmultinom(int*, double*);
extern void rmultinom_unif(int*, int);
extern int rmultinom_unif_fast(int);

extern double rnorm(double*, double*); 
extern double dnorm(double, double*, double*);
extern double pnorm(double, double*, double*, int);
extern double pt(double x, double df, int type);

extern void rdirichlet(double*, double*, int*);
extern double ddirichlet(double*, double*, int*);



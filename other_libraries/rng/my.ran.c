#include "ranlib.h"
#include "cdflib.h"
#include <math.h>

#define PI 3.141592654

double runif() {

	return (double) ranf();
}

double rnorm(double *mu, double *sd) {

	double a;

	a = (double) snorm();

	return *mu + *sd * a;
}

double dnorm(double x, double *mu, double *sd) {

	double a;

	a = exp(-0.5 * (x - *mu) * (x - *mu) / (*sd * *sd)) / (sqrt(2 * PI) * *sd);

	return a;
}

double pnorm(double x, double *mean, double *sd, int type) {
	
	double p, q, z;

	z = (x - *mean) / *sd;

	if(z > 5.0) {	
		p = 0.9999997;
		q = 1.0 - p;
	}
	else if(z < -5.0) {	
		p = 2.86e-7;
		q = 1.0 - p;
	}
	else {	
		cumnor(&z ,&p, &q);
	}
	if(type == 0) return p;

	return q;

}

double pt( double x, double df, int type) {

	double p, q, t, z, df1;
	
	if( df < 10 & x < -3000.0) { p = 0.0; q = 1 - p; }
	else if( df < 10 &  x > 3000.0) { p = 1.0; q = 1 - p; }
	else if ( df >= 10 &  x < -150) { p = 0.0; q = 1 - p; }
	else if ( df >= 10 & x > 150) { p = 1.0; q = 1 - p; }
	else if ( df > 2 ) 
		{
		z = x / sqrt(df) * sqrt( df - 2.0);
		if( df < 10 & z < -3000.0) { p = 0.0; q = 1 - p; }
		else if( df < 10 & z > 3000.0) { p = 1.0; q = 1 - p; }
		else if ( df >= 10 & z < -150) { p = 0.0; q = 1 - p; }
		else if ( df >= 10 & z > 150) { p = 1.0; q = 1 - p; }
		else 
			{
			 z = x; df1 = (double) df;
			 cumt(&z ,&df1, &p, &q);
			}
		}
	else { z = x; df1 = (double) df; cumt(&z ,&df1, &p, &q); }
  	
	if(type == 0) return p;
          return q;
}


double rgamma(double shape, double rate) {
	
	return gengam((float) rate, (float) shape);
}

double dgamma(double x, double a, double b, int flag) {
  
  /* calculates the gamma density */
  /* if flag = 1 it returns the log-density */
  /* if flag = 0 it returns the density */
  /* takes shape(=a) and rate(=b) */

  /* Density is given by f(x)= (b^a / Gamma(a)) x^(a-1) e^(-bx) */

  double d;
  
  d = a * log(b) + (b - 1) * x - b * x - lgamma(a);
  
  if(flag == 0) 
    return exp(d);
  else
    return d;
  	
}

void rmultinom(int *ans, double *probs){

	/* generate a single draw from a Multinomial(probs) distribution */
	int i;
	double unif;
	unif = runif();
	i = 0;
	while(unif > 0) {
		*ans = (i + 1);
		unif -= *(probs + i);
		i += 1;
	}
}

void rdirichlet(double *ans, double *probs, int *l){
	
	/* generate a single draw from a Dirichlet(probs) distribution */
	int i;
	double t = 0.0;
	
	for(i = 0; i < *l; i++) {
		*(ans + i) = rgamma(*(probs + i), 1.0);
		t += *(ans + i);
	}
	for(i = 0; i < *l; i++) *(ans + i) = (*(ans + i) + (0.01 * t)) / (t * (1 + (*l * 0.01)));
	
}

double ddirichlet(double *x, double *pars, int *l) {

	/* returns the log-density of x from a Dirichlet distribution with parameter vector pars */
	int i;
	double t1 = 0.0, log_dens = 0.0;
		
	for(i = 0; i < *l; i++) {
		t1 += *(pars + i);
		log_dens += (*(pars + i) - 1) * log(*(x + i)) - lgamma(*(pars + i));
	}
	
	log_dens += lgamma(t1);
			
	return(log_dens);
}

double lchoose(double n, double x) {

	double a1, a2, a3;
	
	a1 = n + 1;
	a2 = n - x + 1;
	a3 = x + 1;

	return alngam(&a1) - (alngam(&a2) + alngam(&a3));
}


float lchoosef(float n, float x) {

	double a1, a2, a3;
	
	a1 = (double) (n + 1.0);
	a2 = (double) (n - x + 1.0);
	a3 = (double) (x + 1);

	return (float) (alngam(&a1) - (alngam(&a2) + alngam(&a3)));
}

int rmultinom_unif(int n){

	/* generate a single draw from a Multinomial(probs) distribution with uniform probs*/
	int i = 0, ans;
	double unif = runif();
	while(unif > 0) {
		ans = (i + 1);
		unif -= 1.0 / ((double) n);
		i += 1;
	}
	return ans;
}

int rmultinom_unif_fast(int n) {

	/* generate a single draw from a Multinomial(probs) distribution with uniform probs*/
	return (int) floor(((float) n) * runif()) + 1;
}

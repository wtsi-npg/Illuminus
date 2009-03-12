#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_3D_array_int(int *mat, int n, int p, int q) {

	int i, j, k;

	/* printf("\n"); */
	for(i = 0; i < n; i++) {
		for(j = 0; j < p; j++) {
			for(k = 0; k < q; k++) printf("%d ", *(mat + i * p * q + j * q + k));
			printf("\n");
		}
		printf("\n");
	}
	printf("\n");

}

void print_3D_array_double(double *mat, int n, int p, int q) {

	int i, j, k;

	/* printf("\n"); */
	for(i = 0; i < n; i++) {
		for(j = 0; j < p; j++) {
			for(k = 0; k < q; k++) printf("%f ", *(mat + i * p * q + j * q + k));
			printf("\n");
		}
		printf("\n");
	}
	printf("\n");

}

void print_mat_int(int *mat, int n, int p) {

	int i, j;

	/* printf("\n"); */
	for(i = 0; i < n; i++){
		for(j = 0; j < p; j++) printf("%d ", *(mat + i * p + j));
		printf("\n");
	}
	printf("\n");

}

void print_mat_double(double *mat, int n, int p) {

	int i, j;

	/* printf("\n"); */
	for(i = 0; i < n; i++){
		for(j = 0; j < p; j++) printf("%f ", *(mat + i * p + j));
		printf("\n");
	}
	printf("\n");

}

void sequential_mean(double *mu, double *x, double m) {

	/* calculate a mean sequentially */

	*mu = (*x + (m - 1) * *mu) / m;

}

void sequential_variance(double *var, double *mu, double *x, double m) {

	/* calculate the variance sequentially */

	if(m == 1) *var = *x;
	if(m == 2) *var = ((*var - *x) / 2.0) * (*var - *x);
	if(m > 2) *var = (m - 2.0) * *var / (m - 1.0) + m * (*x - *mu) * (*x - *mu) / ((m - 1.0) * (m - 1.0)); 
}

void sequential_sd(double *sd, double *mu, double *x, double m) {

	/* calculate the variance sequentially */

	if(m == 1) *sd = *x;
	if(m == 2) *sd = sqrt(((*sd * *sd - *x) / 2.0) * (*sd * *sd - *x));
	if(m > 2) *sd = sqrt((m - 2.0) * *sd * *sd / (m - 1.0) + m * (*x - *mu) * (*x - *mu) / ((m - 1.0) * (m - 1.0))); 
}


void sequential_mean_vec(double *mu, double *x, double m, int vec_length) {

	int i;

	for(i = 0; i < vec_length; i++) sequential_mean(mu + i, x + i, m);
}

void sequential_variance_vec(double *var, double *mu, double *x, double m, int vec_length) {

	int i;

	for(i = 0; i < vec_length; i++) sequential_variance(var + i, mu + i, x + i, m); 
}

void sequential_sd_vec(double *sd, double *mu, double *x, double m, int vec_length) {

	int i;

	for(i = 0; i < vec_length; i++) sequential_sd(sd + i, mu + i, x + i, m); 
}

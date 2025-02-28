#pragma once
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void print_matrix(const gsl_matrix *m);
void print_vector(const gsl_vector *v);
void rand_tree(gsl_matrix *A,int (*depth)(int n), double(*ro)(int d));
int depth(int n);
double ro(int d);
void compute_eigenvalues(const gsl_matrix *A, gsl_vector *eval, gsl_matrix *evec);
void compute_degree_matrix(const gsl_matrix *A, gsl_matrix *D);





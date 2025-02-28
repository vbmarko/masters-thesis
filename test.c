#include <bits/pthreadtypes.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

void compute_eigenvalues(const gsl_matrix *A, gsl_vector *eval, gsl_matrix *evec) {
    gsl_matrix *A_copy = gsl_matrix_alloc(A->size1, A->size2);
    gsl_matrix_memcpy(A_copy, A);

    gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(A->size1);
    gsl_eigen_symmv(A_copy, eval, evec, workspace);

    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(A_copy);

    for (size_t i = 0; i < eval->size; i++) {
        printf("Eigenvalue %zu: %g\n", i, gsl_vector_get(eval, i));
    }
}


void compute_degree_matrix(const gsl_matrix *A, gsl_matrix *D) {
    size_t n = A->size1;

    // Initialize degree matrix to zeros
    gsl_matrix_set_zero(D);

    for (size_t i = 0; i < n; i++) {
        double degree = 0.0;

        for (size_t j = 0; j < n; j++) {
            degree += gsl_matrix_get(A, i, j);
        }

        gsl_matrix_set(D, i, i, degree); // Set diagonal entry to degree
    }
}

int main(){
	int n = 3;
	gsl_matrix *A = gsl_matrix_alloc(n,n); // Allocate memory for adjacency matrix
    	
	gsl_matrix_set(A, 0,0,0);
    	gsl_matrix_set(A, 1,1,0);
    	gsl_matrix_set(A, 2,2,0);
    	gsl_matrix_set(A, 0,1,1);
   	gsl_matrix_set(A, 1,0,1);
    	gsl_matrix_set(A, 2,1,0);
    	gsl_matrix_set(A, 1,2,0);
    	gsl_matrix_set(A, 0,2,1);
   	 gsl_matrix_set(A, 2,0,1);


	gsl_vector *eval = gsl_vector_alloc(n);
	gsl_matrix *evec = gsl_matrix_alloc(n,n);
	gsl_matrix *L = gsl_matrix_calloc(n,n);
	gsl_matrix *D = gsl_matrix_calloc(n,n);
	 // Compute degree matrix D
    	compute_degree_matrix(A, D);

    	// Compute Laplacian L = D - A
    	gsl_matrix_memcpy(L, D); // L = D
    	gsl_matrix_sub(L, A);    // L = L - A

  printf("Adjacency Matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", (int) gsl_matrix_get(A, i, j));
        }
        printf("\n");
    }


	compute_eigenvalues(L, eval,evec);
	printf("Laplacian:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
	   printf("%d ",(int)gsl_matrix_get(L, i,j));

	}   printf("\n");
    }




	printf("evec:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ",  gsl_matrix_get(evec, i, j));
        }
        printf("\n");
    }

	return 0;
}



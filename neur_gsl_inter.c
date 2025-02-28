#include "nrn_lib/include/oc_ansi.h"
#include "nrn_lib/include/hocdec.h"
#include "nrn_lib/include/hoc.h"
#include <stdio.h>
#include <gsl/gsl_matrix.h>

extern int hoc_valid_stmt(const char* stmt, struct Object* ob);
extern double hoc_ac_;

// Function to extract adjacency matrix based on parent-child relationships
void build_adjacency_matrix(gsl_matrix *A) {
    // Get the number of sections using hoc_valid_stmt
    hoc_valid_stmt("num_sections = 0; forall num_sections = num_sections + 1; hoc_ac_ = num_sections", 0);
    int num_sections_int = (int)hoc_ac_; 
    printf("Number of sections: %d\n", num_sections_int);

    // Clear the adjacency matrix
    gsl_matrix_set_zero(A);

    // Loop through each section and determine its parent
    for (int i = 0; i < num_sections_int; i++) {
        char hoc_parent_command[256];
        snprintf(hoc_parent_command, sizeof(hoc_parent_command),
                 "objref sec_parent; sec_parent = new SectionRef(); hoc_ac_ = sec_parent.has_parent() ? sec_parent.parent.sec_id() : -1");
        hoc_valid_stmt(hoc_parent_command, 0);
        int parent_idx = (int)hoc_ac_;

        if (parent_idx >= 0) {
            gsl_matrix_set(A, i, parent_idx, 1); // Connect child to parent
            gsl_matrix_set(A, parent_idx, i, 1); // Undirected connection
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <path_to_morphology.swc>\n", argv[0]);
        return 1;
    }

    // Initialize NEURON
    hoc_main1(0, (const char **)(argv), 0);

    // Load the provided SWC morphology file
    char hoc_command[512];
    snprintf(hoc_command, sizeof(hoc_command), "load_file(\"%s\")", argv[1]);
    hoc_valid_stmt(hoc_command, 0);

    // Get the number of sections
    hoc_valid_stmt("num_sections = 0; forall num_sections = num_sections + 1; hoc_ac_ = num_sections", 0);
    int num_sections_int = (int)hoc_ac_;

    // Allocate GSL adjacency matrix
    gsl_matrix *A = gsl_matrix_calloc(num_sections_int, num_sections_int);

    // Build adjacency matrix
    build_adjacency_matrix(A);

    // Print adjacency matrix
    printf("Adjacency Matrix:\n");
    for (int i = 0; i < num_sections_int; i++) {
        for (int j = 0; j < num_sections_int; j++) {
            printf("%d ", (int)gsl_matrix_get(A, i, j));
        }
        printf("\n");
    }

    // Free memory
    gsl_matrix_free(A);

    return 0;
}


#include "nrn_lib/include/oc_ansi.h"
#include "nrn_lib/include/hocdec.h"
#include "nrn/src/nrniv/ocjump.h"
#include "nrn_lib/include/ocfunc.h"
//#include <bits/pthreadtypes.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "nrn_lib/include/hocdec.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <cstdarg>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include "./panes/rl_panes.h"
#include "./panes/raylib.h"
#include <string.h>


#include "util.h"
#include "my_nrn_lib.h"
#include "gsl_funs.h"


// extern C functions
extern int hoc_valid_stmt(const char* stmt, struct Object* ob);
extern double hoc_ac_;
// oc environment
static char** env = 0;
void modl_reg(void) { }




int main(int argc, char **argv) {
    hoc_main1(0, (const char **)(argv), (const char **)env);
    NEURON_MEHCANISM *mech = (NEURON_MEHCANISM*) malloc(1*sizeof(NEURON_MEHCANISM));
    init_mechanism(mech);

    print_mech(mech);

    char *template_path = "Cell.hoc";
    ball_n_y_swc("bny.swc", 10, 1, .5, .1, 500, 200, 1);

    char *morph_path = "te2.swc";
    load_file(template_path);
    //load_file(template_path);

    NEURON_CELL *c = (NEURON_CELL*) malloc(sizeof(NEURON_CELL));
    

    double Rm = 3300, v_init = -65;
    init_cell(Rm,v_init,c, "c", morph_path);
    printf("\n");
    print_matrix(c->G);

    printf("\n");


    matrix_to_file(c->G, "G.csv", ',', "%g");
    matrix_to_file(c->sA, "sA.csv", ',', "%g");
    vector_to_file(c->C, "C.csv", ',', "%g");

 


    char *cmd = format_string("%s.C.printf()",c->name);
    hoc_valid_stmt(cmd,0);

    cmd = format_string("%s.print_ri_n_area()",c->name);
    hoc_valid_stmt(cmd, 0);

    gsl_vector *sid = gsl_vector_calloc(c->C->size);

    hoc_to_gsl_vector(sid, c,"seg_id_vec");

    vector_to_file(sid, "sid.csv", ',', "%g");




    return 0;





    double dt = 0.025;     /* time step in ms */
    double tstop = 5.0;  /* simulation end time in ms */
    int nsteps = (int)(tstop / dt);
    return 0;
}


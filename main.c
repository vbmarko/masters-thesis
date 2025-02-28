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
    char *morph_path = "bny.swc";
    load_file(template_path);
    //load_file(template_path);
    //ball_n_y_swc("bny.swc", 8, 1, .2, .1, 500, 250, 5);

    NEURON_CELL *c = (NEURON_CELL*) malloc(sizeof(NEURON_CELL));
    

    double Rm = 33, v_init = 0;
    init_cell(c, "c", morph_path);
    print_matrix(c->G);
    return 0;





    double dt = 0.025;     /* time step in ms */
    double tstop = 5.0;  /* simulation end time in ms */
    int nsteps = (int)(tstop / dt);
    return 0;
}


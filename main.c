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
    char *morph_path = "bny.swc";
    load_file(template_path);
    //load_file(template_path);
    //ball_n_y_swc("bny.swc", 8, 1, .2, .1, 500, 250, 5);

    NEURON_CELL *c = (NEURON_CELL*) malloc(sizeof(NEURON_CELL));
    

    double Rm = 3300, v_init = -65;
    init_cell(Rm,v_init,c, "c", morph_path);
    printf("\n");
    //print_matrix(c->G);

    char *cmd = format_string("%s.matrix_to_svc(%s,%s)",c->name,"c.G","\"G.svc\"");
    hoc_valid_stmt(cmd,0);

    cmd = format_string("%s.vector_to_svc(%s,%s)",c->name,"c.C","\"C.svc\"");
    hoc_valid_stmt(cmd,0);

    printf("\n\n %f",c->D->data[0]);
    return 0;





    double dt = 0.025;     /* time step in ms */
    double tstop = 5.0;  /* simulation end time in ms */
    int nsteps = (int)(tstop / dt);
    return 0;
}


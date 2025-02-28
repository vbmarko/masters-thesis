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


    char *template_path = "Cell.hoc";
    char *morph_path = "c91662.swc";
    //load_file(template_path);
    load_file(template_path);
    NEURON_CELL *c = (NEURON_CELL*) malloc(sizeof(NEURON_CELL));
    init_cell(c, (char*)"c", morph_path);
  


    printf("now lets try to insert a chan %s  at section %s\n",mech->chan[2],c->sec[0].name);
    insert_chan(mech->chan[2]);
    char *cmd = format_string("access %s",c->sec[0].name);
    hoc_valid_stmt(cmd,0);
    hoc_valid_stmt("psection()", 0);


    double s = 0.5;
    printf("now lets try to insert a synapse %s  at section %s at positon %f\n",mech->syn[0],c->sec[0].name,s);
    insert_syn_at_sec(mech->syn[0],"IClamp5", 0, c,  s);

    printf("the cell %s now have %d synapses",c->name,c->n_syn);
    printf("itsa a %s named  %s and its inserted at %f in sec %d",c->syn[0].type,c->syn[0].name,c->syn[0].s,c->syn[0].sec_id );
    
	


 
     return 0; 

    double dt = 0.025;     /* time step in ms */
    double tstop = 5.0;  /* simulation end time in ms */
    int nsteps = (int)(tstop / dt);
    
    /* Allocate a buffer to store the soma voltage at each time step */
    double *voltage_buffer = (double *)malloc(nsteps * sizeof(double));
    if (!voltage_buffer) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

   
    	for (int i = 0; i < nsteps; i++) {
	hoc_valid_stmt("fadvance()",0);
        
        if (hoc_valid_stmt("hoc_ac_ = c.soma.v(0.5)", 0) < 0) {
            fprintf(stderr, "Error evaluating soma.v(0.5) at step %d\n", i);
            voltage_buffer[i] = 0.0;
        } else {
            voltage_buffer[i] = hoc_ac_;
        }
    }
    
    /* Optionally, print out the results */
    for (int i = 0; i < nsteps; i++) {
        printf("t = %.3f ms, V_soma = %.3f mV\n", i * dt, voltage_buffer[i]);
    }

    return 0;
}


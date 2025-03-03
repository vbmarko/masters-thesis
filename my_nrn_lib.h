#pragma once
#include <gsl/gsl_vector_double.h>
#include <gsl//gsl_matrix.h>
typedef struct NEURON_SEGMENT NEURON_SEGMENT;
typedef struct NEURON_SECTION NEURON_SECTION;
typedef struct NEURON_CELL NEURON_CELL;
typedef struct NEURON_MECHANISM NEURON_MEHCANISM;
typedef struct NEURON_SYNAPSE NEURON_SYNAPSE;


struct NEURON_SEGMENT{
	int sec_id;
	double s; //relative position in section
	double v;
	double *g;
	double **mech_state;

};

struct NEURON_SECTION{
	char *name;
	int id, p_id, n_seg, n3d;
	int *n_mech_state;
	double L, Ra, diam;
	double* arc3d;
	double* x;
	double* y;
	double* vx;
	double* vy;


};

struct NEURON_CELL{
	char* name;
	int n_sec, n_syn;
	NEURON_SECTION *sec;
	NEURON_SYNAPSE *syn;
	double *chi_x;
	double *chi_y;
	gsl_matrix *G;
	gsl_matrix *A;
	gsl_matrix *L;
	gsl_matrix *Lw;
	gsl_matrix *D;
	gsl_vector *C;
	gsl_vector *d;
	gsl_vector *eval;
	gsl_matrix *evec;
	
};

struct NEURON_MECHANISM{
	int n_chan, n_syn, n_acell, n_misc;
	char **chan; 
	char **syn;
	char **acell;
	char **misc;

};


struct NEURON_SYNAPSE{
	int sec_id;
	char *type;
	char* name;
	double del;
	double dur;
	double amp;
	double gmax;
	double e;
	double s;
};


void free_cell(NEURON_CELL *c);
void init_mechanism(NEURON_MEHCANISM *mech);
void init_cell(double Rm, double v0, NEURON_CELL *c, char *cell_name, char* morph_path);
void insert_chan_at_sec(char *chan,int i_sec,NEURON_CELL *c);
void insert_syn_at_sec(char *syn,char *syn_name,int i_sec,NEURON_CELL *c,double s);
void insert_chan(char *chan);
int paren_idx(char *syn, NEURON_CELL *c);
void print_sec(NEURON_SECTION *sec);
void update_sec_embedding(int fix_idx,NEURON_SECTION sec, double f_0, double f_1,double k);
void ball_n_y_swc(char *fpath,double soma_diam,double diam_i,double diam_b,double diam_t,double L_stick, double L_branch,double res);
void print_mech(NEURON_MEHCANISM *m);
void set_cell_properties(NEURON_CELL *c,double Rm,double v_init);

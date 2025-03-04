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
#include <math.h>
#include "util.h"
#include "my_nrn_lib.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include "gsl_funs.h"



extern double hoc_ac_;
extern int hoc_valid_stmt(const char* stmt, struct Object* ob);


int matrix_ncol(NEURON_CELL *c,char *m_name){
		char *cmd = format_string("hoc_ac_ = %s.%s.nrow()",c->name,m_name);
		hoc_valid_stmt(cmd, 0);
		return (int) hoc_ac_;
}

int matrix_nrow(NEURON_CELL *c,char *m_name){
		char *cmd = format_string("hoc_ac_ = %s.%s.nrow()",c->name,m_name);
		hoc_valid_stmt(cmd, 0);
		return (int) hoc_ac_;
}


void hoc_to_gsl_matrix(gsl_matrix* m,NEURON_CELL *c, char *m_name){
	for(int i = 0;i<m->size1;++i){
		for(int j = 0;j<m->size2;++j){
			char *cmd = format_string("hoc_ac_ = %s.%s.x[%d][%d]",c->name,m_name,i,j);
			hoc_valid_stmt(cmd, 0);
			m->data[i*m->tda + j] = hoc_ac_;
		}
	}

}

void hoc_to_gsl_vector(gsl_vector* v,NEURON_CELL *c, char *v_name){
	for(int i = 0;i<v->size;++i){
			char *cmd = format_string("hoc_ac_ = %s.%s.x[%d]",c->name,v_name,i);
			hoc_valid_stmt(cmd, 0);
			gsl_vector_set(v,i,hoc_ac_);
		}
	}


void free_cell(NEURON_CELL *c) {
    if (c != NULL) {
        // Free each section's dynamically allocated members
        if (c->sec != NULL) {
            for (int i = 0; i < c->n_sec; i++) {
                // Free the section name allocated in init_cell
                if (c->sec[i].name != NULL) {
                    free(c->sec[i].name);
                    c->sec[i].name = NULL;
                }
                // If you later allocate other members in NEURON_SECTION
                // (for example, c->sec[i].n_mech_state) then free them here.
            }
            // Free the array of sections itself
            free(c->sec);
            c->sec = NULL;
        }
        // Free the NEURON_CELL structure
        free(c);
    }
}


int paren_idx(char *sec,NEURON_CELL *c){
	char *line = (char*)calloc(256,sizeof(char));
	char *cmd = "objref sr";
	hoc_valid_stmt(cmd,0);
	//cmd = format_string("%s.indexSections()",c->name);
	hoc_valid_stmt(cmd,0);
	cmd = format_string("access %s",sec);
	hoc_valid_stmt(cmd,0);
	cmd = "sr = new SectionRef()";
	hoc_valid_stmt(cmd,0);
	cmd = "if(sr.has_trueparent){sr.trueparent {print secname()}}";
	capture_out((char*)line,cmd);
	size_t len = strlen(line);
	line[len-1] = '\0';
	for(int i = 0;i<c->n_sec;++i){
		if(strcmp((char*)line, c->sec[i].name) == 0){
			printf("the sec %s has parent %s \n",c->sec[i].name,line);
			return i;
		}
	}
	return -1;
}



void init_cell(double Rm, double v0, NEURON_CELL *c, char *cell_name, char* morph_path){
	c->n_syn = 0;
	char *cmd = format_string("objref %s",cell_name);
	hoc_valid_stmt(cmd,0);
	cmd = format_string("%s = new Cell()",cell_name);
	hoc_valid_stmt(cmd,0);
	cmd = format_string("%s.AddMorph(\"%s\")",cell_name,morph_path);
	hoc_valid_stmt(cmd,0);

	
	cmd = format_string("%s.delet_axon()",cell_name);

	cmd = format_string("%s.SetCellProperties(%g,%g)",cell_name,Rm,v0);
	hoc_valid_stmt(cmd,0);

	cmd = format_string("%s.geom_nseg()",cell_name);
	hoc_valid_stmt(cmd,0);

	cmd = format_string("print %s.nSecAll",cell_name);
	hoc_valid_stmt(cmd,0);
	cmd = format_string("hoc_ac_ = %s.nSecAll",cell_name);
	hoc_valid_stmt(cmd,0);
	c->n_sec = hoc_ac_;
	c->name = cell_name;
	c->sec = (NEURON_SECTION*) malloc(sizeof(NEURON_SECTION)*c->n_sec);



	cmd = format_string("%s.init_G_n_C()",cell_name);
	hoc_valid_stmt(cmd,0);
	cmd = format_string("%s.init_D()",cell_name);
	hoc_valid_stmt(cmd,0);
	int n_r = matrix_nrow(c, "G"), n_c = matrix_ncol(c, "G");
	gsl_matrix *G = gsl_matrix_calloc(n_c,n_r);
	gsl_matrix *D = gsl_matrix_calloc(n_c,n_r);
	hoc_to_gsl_matrix(G,c, "G");
	hoc_to_gsl_matrix(D,c, "D");


	gsl_matrix *As = gsl_matrix_calloc(c->n_sec, c->n_sec);
	//cmd = format_string("%s.init_sec_adj()",cell_name);
	//hoc_valid_stmt(cmd,0);

	
	hoc_to_gsl_matrix(As,c, "sec_adj_mat");
	c->sA = As;


	gsl_matrix *A = gsl_matrix_calloc(n_c,n_r);
	binarize_matrix(A,G,0);
	
	gsl_matrix *L = gsl_matrix_calloc(n_c,n_r);
	gsl_matrix *Lw = gsl_matrix_calloc(n_c,n_r);
	compute_laplacian_matrix(L,A);
	compute_laplacian_matrix(Lw,G);

	gsl_matrix *evec = gsl_matrix_calloc(n_c,n_r);
	gsl_vector *eval = gsl_vector_calloc(n_c);
	compute_eigenvalues(eval,evec,A);

	gsl_vector *d = gsl_vector_calloc(n_c);
	compute_degree_vector(d,A);

	gsl_vector* C = gsl_vector_calloc(n_c);
	hoc_to_gsl_vector(C, c, "C");

	c->G = G;
	c->C = C;
	c->D = D;
	c->L = L;
	c->Lw = Lw;
	c->d = d;
	c->A = A;
	c->evec = evec;
	c->eval = eval;

	

	//Awfull jank below
	char *line = (char *)calloc(256,sizeof(char));
        hoc_valid_stmt("objref sr",0);
	load_file("print_ith_sec.hoc");
	for(int i = 0;i<c->n_sec;++i){
		cmd = format_string("print_ith_sec(%d,%s,%s)",i,cell_name,"sr");
		capture_out(line,cmd);
		size_t len = strlen(line);
		c->sec[i].name = (char *) calloc(len+1,sizeof(char));
		strcpy(c->sec[i].name, line);
		c->sec[i].name[len-1] = '\0';
		c->sec[i].id = i;
				
		cmd = format_string("access %s",c->sec[i].name);
		hoc_valid_stmt(cmd, 0);
		//hete we add L, Ra and diam
		cmd = "hoc_ac_ = L";
		hoc_valid_stmt(cmd, 0);
		
		c->sec[i].L = hoc_ac_;
		cmd = "hoc_ac_ = Ra";
		hoc_valid_stmt(cmd, 0);
		c->sec[i].Ra = hoc_ac_;
		cmd = "hoc_ac_ = diam";
		hoc_valid_stmt(cmd, 0);
		c->sec[i].diam = hoc_ac_;

			
		//here we add 3dpnts
	
		cmd = "hoc_ac_ = n3d()";
		hoc_valid_stmt(cmd, 0);
		int n = hoc_ac_;
		c->sec[i].n3d = n;
		c->sec[i].arc3d = (double*) calloc(n,sizeof(double));
		c->sec[i].x = (double*) calloc(n,sizeof(double));
		c->sec[i].y = (double*) calloc(n,sizeof(double));
		for(int j = 0;j<c->sec[i].n3d;++j){
			cmd =format_string("hoc_ac_ = arc3d(%d)",j);
			hoc_valid_stmt(cmd, 0);
			c->sec[i].arc3d[j] = hoc_ac_;
			c->sec[i].x[j] = 1.0*rand()/RAND_MAX;
			c->sec[i].y[j] = 1.0*rand()/RAND_MAX;

		}
	


	}
	for(int i = 0;i<c->n_sec;++i){
			int idi = paren_idx(c->sec[i].name,c);
			c->sec[i].p_id = idi;
			printf("hello from init neuron the secs idx is %d the parents idx is %d\n",c->sec[i].id,idi);
		}
	free(line);
	printf("Hello fren the cell is initalized please free it with free_cell after youre done with it!");
}


void insert_chan_at_sec(char *chan,int i_sec,NEURON_CELL *c){
	char *cmd = format_string("access %s",c->sec[i_sec].name);
	hoc_valid_stmt(cmd, 0);
	cmd = format_string("insert %s",chan);
	hoc_valid_stmt(cmd, 0);
}

void insert_syn_at_sec(char *syn,char *syn_name,int i_sec,NEURON_CELL *c,double s){
	char *cmd = format_string("access %s",c->sec[i_sec].name);
	hoc_valid_stmt(cmd, 0);
	cmd = format_string("objref %s",syn_name);
	hoc_valid_stmt(cmd, 0);
	cmd = format_string("%s %s = new %s(0.5)",c->sec[i_sec].name,syn_name,syn);
	hoc_valid_stmt(cmd, 0);
	hoc_valid_stmt(cmd,0);
	c->syn = (NEURON_SYNAPSE*)realloc(c->syn, c->n_syn + 1);
	c->syn[c->n_syn] = (NEURON_SYNAPSE) {.sec_id = i_sec,.type = syn,.name = syn_name,.s = s};
	c->n_syn += 1;
	
}




void insert_chan(char *chan){
	char *cmd = format_string("forall {insert %s}",chan);
	hoc_valid_stmt(cmd, 0);
}


void init_mechanism(NEURON_MEHCANISM *mech){
	load_file("safe_mech_count.hoc");
	hoc_valid_stmt("strdef mname", 0);
	char *line = (char*) malloc(256);
	char *cmd = "objref mt";
	hoc_valid_stmt(cmd, 0);
	cmd = format_string("mt = new MechanismType(%d)",0);
	hoc_valid_stmt(cmd, 0);
	cmd = format_string("hoc_ac_ = mt.count()");
	hoc_valid_stmt(cmd, 0);
	mech->n_chan = hoc_ac_;
	mech->chan = (char**) malloc(mech->n_chan*sizeof(char*));

	for(int i = 0;i<mech->n_chan;++i){
		cmd = format_string("mt.select(%d)",i);
		hoc_valid_stmt(cmd, 0);
		cmd = format_string("mt.selected(mname)",i);
		hoc_valid_stmt(cmd, 0);
		cmd = "print mname";
		capture_out(line,cmd);
		size_t len = strlen(line);
		mech->chan[i] = (char *) calloc(len+1,sizeof(char));
		strcpy(mech->chan[i], line);
		mech->chan[i][len-1] = '\0';

	}	

	cmd = format_string("mt = new MechanismType(%d)",1);
	hoc_valid_stmt(cmd, 0);
	cmd = format_string("hoc_ac_ = mt.count()");
	hoc_valid_stmt(cmd, 0);	
	mech->n_syn = hoc_ac_;
	mech->syn = (char**) malloc(mech->n_syn*sizeof(char*));
	for(int i = 0;i<mech->n_syn;++i){
		cmd = format_string("mt.select(%d)",i);
		hoc_valid_stmt(cmd, 0);
		cmd = format_string("mt.selected(mname)",i);
		hoc_valid_stmt(cmd, 0);		
		cmd = "print mname";
		capture_out(line,cmd);
		size_t len = strlen(line);
		mech->syn[i] = (char *) calloc(len+1,sizeof(char));
		strcpy(mech->syn[i], line);
		mech->syn[i][len-1] = '\0';



	}


/*
	cmd = format_string("mt = new MechanismType(%d)",2);
	hoc_valid_stmt(cmd, 0);
	cmd = format_string("hoc_ac_ = mt.count()");
	hoc_valid_stmt(cmd, 0);	
	mech->n_acell = hoc_ac_;
	mech->acell = (char**) malloc(mech->n_acell*sizeof(char*));
	for(int i = 0;i<mech->n_acell;++i){
		cmd = format_string("mt.select(%d)\nmt.selected(mname)",i);
		hoc_valid_stmt(cmd, 0);
		cmd = "print mname";
		capture_out(line,cmd);
		size_t len = strlen(line);
		mech->acell[i] = (char *) calloc(len+1,sizeof(char));
		strcpy(mech->acell[i], line);


	}	
	cmd = format_string("hoc_ac_ = safe_mech_count(%d,%s)",3,"mt");
	hoc_valid_stmt(cmd, 0);
	mech->n_misc = hoc_ac_;
	mech->misc = (char**) malloc(mech->n_misc*sizeof(char*));
	

	for(int i = 0;i<mech->n_misc;++i){
		cmd = format_string("mt.select(%d)\nmt.selected(mname)",i);
		hoc_valid_stmt(cmd, 0);
		cmd = "print mname";
		capture_out(line,cmd);
		size_t len = strlen(line);
		mech->misc[i] = (char *) calloc(len+1,sizeof(char));
		strcpy(mech->misc[i], line);
	}
*/	
	free(line);
}



void print_mech(NEURON_MEHCANISM *m){
	printf("Available channels:\n");
	for(int i = 0;i<m->n_chan;++i){
		printf("%s \n",m->chan[i]);
	}	
	printf("Available synapses:\n");
	for(int i = 0;i<m->n_syn;++i){
		printf("%s \n",m->syn[i]);
	}

	
}



void print_sec(NEURON_SECTION *sec){
	printf("\n\nHello I am a section named %s my index is %d my parents index is %d\n",sec->name,sec->id,sec->p_id);
	printf("My biophisical properties are: \n");
	printf("L = %.2f, Ra = %.2f, diam = %.2f\n",sec->L,sec->Ra,sec->diam);
	printf("The physiologist measured me %d points along my extent: \n",sec->n3d);
	for(int i = 0;i<sec->n3d;++i){
		printf("	%d at %.2f\n",i,sec->arc3d[i]);
	}
}



void update_sec_embedding(int fix_idx,NEURON_SECTION sec, double f_0, double f_1,double k){
	int n = sec.n3d;
	double err = 0;
	double vx[n], vy[n];
	for(int i = 0;i<n;++i){
		double x = sec.x[i], y = sec.y[i], s = sec.arc3d[i];
		if(i == fix_idx){
			vx[i] = 0;
			vy[i] = 0;
		}else if(i == 0){
	
			double lf = fabs(s - sec.arc3d[i+1]);
			double dxf = x - sec.x[i+1], dyf = y- sec.y[i+1];
			double df = sqrt(dyf*dyf + dxf*dxf);
			vx[i] = (k)*(lf-df)*dxf  + f_0;
			vy[i] = (k)*(lf-df)*dyf  + f_0;
			err += fabs(lf-df);


			
		}else if (i == n-1) {

			double lb = fabs(s - sec.arc3d[i-1]);
			double dxb = x -sec.x[i-1], dyb = y - sec.y[i-1];
			double db = sqrt(dxb*dxb + dyb*dyb);
			vx[i] = (k)*(lb-db)*dxb  + f_1;
			vy[i] = (lb-db)*db*dyb  + f_1;
			err += fabs(lb-db);


			
		}else{
			
			double lb = fabs(s - sec.arc3d[i-1]), lf = fabs(s - sec.arc3d[i+1]);
			double dxb = x -sec.x[i-1], dyb = y - sec.y[i-1], dxf = x - sec.x[i+1], dyf = y- sec.y[i+1];
			double db = sqrt(dxb*dxb + dyb*dyb), df = sqrt(dyf*dyf + dxf*dxf);
			double x_dot = (k)*(lb-db)*dxb  + (k)*(lf-df)*dxf;
			printf(" the x dot %.2f \n",x_dot);
			vx[i] = x_dot;
			vy[i] = (k)*(lb-db)*dyb  + (k)*(lf-df)*dyf;
			err += fabs(lb-db) + fabs(lf-df);

			
		}
	}
	printf("\n----------------------\nhello from embedding of sec %s the err currently is %.2f\n",sec.name,err);
	printf("and the positons\n");
	for(int i = 0;i<sec.n3d;++i){
		sec.x[i] += vx[i];
		sec.y[i] += vy[i];
		printf("(x,y) = (%.2f,%.2f) and real arclen s = %.2f and apparent arclen\n",sec.x[i],sec.y[i],sec.arc3d[i]);
	}

}

void ball_n_y_swc(char *fpath,double soma_diam,double diam_i,double diam_b,double diam_t,double L_stick, double L_branch,double res){
	int N_s = (int) 1, N_i = (int) L_stick*res, N_y = (int) L_branch*res;
	char* head = "#swc with an apical dendrite with one branching of eqaul lengths\n";
        char *templ = "%d %d %f %f %f %f %d\n";  //"1 1 0.0 0.0 0.0 8.8677 -1" = id type x y z radius p_id
	file f = toggle_capture_out_on(fpath);

	printf("%s",head);
	double x = 0, y = 0, z = 0, diam = soma_diam;
	int id = 0, p_id = -1, type = 1;
	//soma
	printf(templ,id,1,x,y,z,diam,p_id);
	type = 4;
	double t = 0;
	for(int i = 0;i<=N_i;++i){
		t = 1.0*i/(N_i);
		id += 1;
		p_id += 1;
		y += L_stick*(1)/(N_i);
		diam = diam_i*(1-t) + diam_b*t;  
		printf(templ,id,type,x,y,z,diam,p_id);
	}
	double normr = sqrt(2);
	double xl = 0, yl = y;
	double b_id = id;
	type = 4;
	for(int i = 0;i<N_y;++i){
		t = 1.0*(i)/(N_y);
		id += 1;

		p_id += 1;
		xl += (L_branch/normr)*(1)/(N_y);
		yl += (L_branch/normr)*(1)/(N_y);
		diam = diam_b*(1-t) + diam_t*t;  
		printf(templ,id,type,xl,yl,z,diam,p_id);

	}	
	type = 4;
	xl = 0, yl =y;
	for(int i = 0;i<N_y;++i){
		if(i == 0){
			p_id = b_id;
		}else if (i == 1) {
			p_id = id;
		
		}
		t = 1.0*(i)/(N_y);
		id += 1;

		xl -= (L_branch/normr)*(1)/(N_y);
		yl += (L_branch/normr)*(1)/(N_y);
		diam = diam_b*(1-t) + diam_t*t;  
		printf(templ,id,type,xl,yl,z,diam,p_id);

		if(i != 0){p_id += 1;};
	}

	toggle_capture_out_off(&f);

	
}



void set_cell_properties(NEURON_CELL *c,double Rm,double v_init){
	char *cmd = format_string("c.SetCellProperties(%f,%f)",Rm,v_init);
	hoc_valid_stmt(cmd, 0);
}






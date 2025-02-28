#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "raylib.h"
#include <math.h>
#include <omp.h>
#include <stdarg.h>

double id( double x){return x;}


//ENTERING
//STRUCTLAND
//


//computing structs
typedef struct net net;
typedef struct compartment compartment;
typedef struct channel channel;
typedef struct gating_params gating_params ;

//preprocesing sphagetti 
typedef struct ptr_w_len ptr_w_len ;
typedef struct segment segment;
typedef struct point point;

//preprocesing sphagetti definitions 
struct ptr_w_len{
 	void *ptr;
	int len;
};
struct point{
	int id;
	int type;
	double x;
	double y;
	double z;
	double r;
	int parent_id;
	int *child_ids;
	int n_childs;
};
struct segment{
	double l;
	double r_i;
	double r_t;
	double r_m;
	double x;
	double y;
	double z;
	double dx;
	double dy;
	double dz;
	segment *p;
	segment **c;
	int n_c;
	int type;
	int id;
	int p_id;
	int *c_id;
	double *u;
	double *u_;
	double *q;
	double g_p;
	double *g_c;
	double n_u;
	double v;
	double v_;
	double S;

};


//computing struct definitions
struct compartment{
	int id;
	int n_ne;
	int *id_ne;
	double *W;
	int n_W;
	double *g_ne;
	double *g_max;
	double *g_s;
	double S;
	double v;
	double temp_v;
	double tau;
	double d;
	int type;
	int n_chans;
	channel *chans;

};struct gating_params{
    double V_m;    // Midpoint voltage (mV)
    double k;      // Slope factor (mV)
    double scale;  // Scaling factor
};


struct channel{
	int n_m;
	double *m;
	double *tau;
	int *exponents;
	gating_params *m_params;
	double (*m_inf)(double v);
	double reversal;

};struct net{
	int n;
	compartment *cmp;
	double dt;
	void (*I_m)(compartment *cmp,double dt);
	void (*I_l)(compartment *cmp, double dt);
	void (*upd_cmp)(compartment *cmp);
	void (*upd_chan)(compartment *cmp,channel *chan,double dt);
	void (*upd)(net *net);

};


//freeland 
//of structs // grace be to chat gpt

// Free functions for specific types
void free_point(point *pt) { free(pt->child_ids); }
void free_segment(segment *seg) {
    free(seg->c);
    free(seg->c_id);
    free(seg->g_c);
    free(seg->u);
    free(seg->u_);
    free(seg->q);
}
void free_channel(channel *channel) {
    free(channel->m);
    free(channel->tau);
    free(channel->exponents);
    free(channel->m_params);
}
void free_compartment(compartment *cmp) {
    free(cmp->id_ne);
    free(cmp->g_ne);
    free(cmp->g_max);
    free(cmp->g_s);
    for (int i = 0; i < cmp->n_chans; i++) {
        free_channel(&cmp->chans[i]);
    }
    free(cmp->chans);
}
void free_net(net *net) {
    for (int i = 0; i < net->n; i++) {
        free_compartment(&net->cmp[i]);
    }
    free(net->cmp);
}

// Macro to define cleanup types
#define CLEANUP_POINT 1
#define CLEANUP_SEGMENT 2
#define CLEANUP_COMPARTMENT 3
#define CLEANUP_NET 4

// Variadic cleanup function
void cleanup_resources(int count, ...) {
    va_list args;
    va_start(args, count);

    for (int i = 0; i < count; i++) {
        int type = va_arg(args, int);
        void *ptr = va_arg(args, void *);

        switch (type) {
            case CLEANUP_POINT:
                free_point((point *)ptr);
                free(ptr);
                break;
            case CLEANUP_SEGMENT:
                free_segment((segment *)ptr);
                free(ptr);
                break;
            case CLEANUP_COMPARTMENT:
                free_compartment((compartment *)ptr);
                free(ptr);
                break;
            case CLEANUP_NET:
                free_net((net *)ptr);
                free(ptr);
                break;
            default:
                fprintf(stderr, "Unknown type for cleanup\n");
        }
    }

    va_end(args);
}


//
//
//
//  ENTERING
//  PRINTLAND

void print_cmp(compartment *cmp){
	printf("id = %d\n",cmp->id);
	printf("n_ne = %d\n",cmp->n_ne);
	printf("S = %f\n",cmp->S);
	printf("v = %f\n",cmp->v);
	printf("----------------------------\n");
}

void print_net(net *net,int n){
	for(int i = 0;i<n;++i){
		compartment *cmp = (net->cmp + i);
		printf("id = %d\n",cmp->id);
		printf("type = %d\n",cmp->type);
		printf("n_ne = %d\n",cmp->n_ne);
		printf("id_ne = [");
		for(int k = 0;k<cmp->n_ne;++k){
			printf("%d,",cmp->id_ne[k]);
		}	

		printf("]\n");
		printf("d = %f\n",cmp->d);
		printf("g_ne = [");
		for(int k = 0;k<cmp->n_ne;++k){
			printf("%f,",cmp->g_ne[k]);
		}
		printf("]\n");
		printf("S = %f\n",cmp->S);
		printf("v = %f\n",cmp->v);
		printf("temp_v=%f\n",cmp->temp_v);
		printf("n_chans = %d\n",cmp->n_chans);
		printf("g_max = [");
		for(int k = 0;k<cmp->n_chans;++k){
			printf("%f,",cmp->g_max[k]);
		}
		printf("]\n-------------------------\n");
	}
}


//
//
//
// ENTERING
// PREPROCESINGLAND


//here we compute the conductances between segments and add them to the segs structs
void add_gs(segment *segs,int n,double r_l){
	for(int i =0;i<n;++i){
		segment *seg = segs + i;
		segment *paren = seg->p;
		double R_i = (r_l * seg->l ) / (2*PI*seg->r_m*seg->r_i);
		double R_p = (r_l * paren->l ) / (2*PI*paren->r_m*paren->r_t);
		seg->g_p = 1/((seg->S)*(R_i + R_p));
		for(int k = 0;k < seg->n_c;++k){
			segment *child = *(seg->c+k);
			double R_i = (r_l * seg->l / 2) / (PI*seg->r_m*seg->r_t);
			double R_p = (r_l * child->l /2) / (PI*child->r_m*child->r_i);
			seg->g_c[k] = 1/((seg->S)*(R_i + R_p));

		}

		
		
	}

}




// here we initialize a segments from points 
segment *new_segments(point *pts, int n){
	segment *segs = malloc(n*sizeof(segment));
	for(int i = 0;i<n;++i){
		point *pt_i = (pts+i);
			
		segment *seg = segs + pt_i->id;
		seg->x = pt_i->x, seg->y = pt_i->y, seg->z = pt_i->z; 
		seg->c_id = pt_i->child_ids;	
		if(pt_i->parent_id < 0){
			seg->type = 0; seg->r_i = pt_i->r; seg->r_t = pt_i->r; seg->l = 0;
			seg->r_m = (seg->r_i+seg->r_t)/2 ; 
			seg->p = seg;
			seg->c = (segment**)(NULL);
			seg->n_c = 0;
			seg->id = pt_i->id;
			seg->v = 10.0 * rand() / RAND_MAX;
			seg->v_ = 0.0;
			seg->S = 4*PI*seg->r_i*seg->r_i;
			seg->l = seg->r_i;
		
		}else{
			point *paren = pts + pt_i->parent_id;
			segment *paren_seg = segs + pt_i->parent_id;
			seg->p_id = pt_i->parent_id;
				
			double dx = paren->x - pt_i->x,dy = paren->y - pt_i->y,dz = paren->z - pt_i->z;


			seg->type = 1; seg->r_i = paren->r; seg->r_t = pt_i->r; seg->l = sqrt(dx*dx+dy*dy+dz*dz);
			seg->r_m = (seg->r_i +seg->r_t) / 2;

			seg->p = segs + pt_i->parent_id;
			if(seg->p->id == 0){
				seg->r_i = seg->r_t;
			}
			seg->S = PI*seg->l*(seg->r_i+seg->r_t);
			seg->c = (segment**)(NULL);
			seg->n_c = 0;
			seg->dx = dx;seg->dy = dy;seg->dz=dz;

			seg->id = pt_i->id;

			paren_seg->n_c = paren_seg->n_c + 1;
			paren_seg->c = realloc(paren_seg->c, sizeof(segment*)*paren_seg->n_c);
			paren_seg->g_c = realloc(paren_seg->g_c, sizeof(double*)*paren_seg->n_c);
			paren_seg->c[paren_seg->n_c-1] = seg;
			seg->v = 10.0 * rand() / RAND_MAX;
			seg->v_ = 0.0;
		}

	}
	return segs;

}





char **break_at_delim(char *str,char delim){
	int toto_sz = 1;
	int num_lines = 1;
	char **lines = malloc(num_lines*sizeof(char **));
	while(*str != '\0'){
		
		int line_sz = 1;
		char *line = (char *)malloc(line_sz);
		while(*str != delim){
			line[line_sz-1] = *str;
			++line_sz;
			line = realloc(line, line_sz);
			++str;
			

		}
		line[line_sz-1] = '\0';
		lines[num_lines-1] = line;
		//free(line);
		++num_lines;
		lines = realloc(lines, sizeof(char *)*num_lines);
		++str;
	}
lines[num_lines-1] = "\0"; 
return lines;

}


void init_segs(char **lines,int n,point *memory,double scale){
	point segs[n];
	for(int i=0;i<n;++i){
		char **sublines = break_at_delim(*(lines+i),' ');
		printf("%s\n",*(lines+i));
		printf("n childs = %d \n",memory[i].n_childs);
		memory[i] = (point){.id = (int)atoi(*(sublines+0)) - 1,
				    .type = atoi(*(sublines+1)),
				    .x = scale*atof(*(sublines+2)),
				    .y = scale*atof(*(sublines+3)),
				    .z = scale*atof(*(sublines+4)),
				    .r = scale*atof(*(sublines+5)),
				    .parent_id = atoi(*(sublines+6)) - 1
				     };

		int parent_id = memory[i].parent_id, my_id = memory[i].id;
		if(parent_id < 0){
			continue;
		}else{
			memory[parent_id].n_childs = memory[parent_id].n_childs + 1;
			memory[parent_id].child_ids = realloc(memory[parent_id].child_ids, memory[parent_id].n_childs);
			memory[parent_id].child_ids[memory[parent_id].n_childs - 1] = my_id;
		}

				    
	}
}



void DrawSegments(segment *segs,compartment *cmps, int n,double min_v,double max_v,int cur_i){
	for(int i=0;i<n;++i){
		segment *seg = segs + i ;
		compartment *cmp = cmps + i; 
		double clipd_v = cmp->v*(cmp->v < max_v && cmp->v > min_v) + max_v*(cmp->v > max_v) + min_v*(cmp->v < min_v);
		double ratio = (clipd_v-min_v)/(max_v-min_v);
		Vector3 seg_pos = (Vector3){seg->x ,seg->y ,seg->z };
		Color c = (Color){.r = (int)(255*ratio), .g = 0,.b = (int)255*(1-ratio),.a = 255};
	
		if(seg->id == 0){
			DrawSphere(seg_pos, seg->r_i, c);	
			if(i == cur_i){
			c = BLACK;
			DrawSphereWires(seg_pos,1.1 * seg->r_i,5,5,GREEN);				
		}


		} else{
			segment *paren = seg->p;
			Vector3 paren_pos = (Vector3){paren->x ,paren->y ,paren->z };
			DrawCylinderEx(seg_pos, paren_pos, seg->r_t ,paren->r_i , 5, c);	
			if(i == cur_i){
			c = BLACK;
			DrawCylinderWiresEx(seg_pos,paren_pos,1.1 * seg->r_t, 1.1 * paren->r_i ,5,GREEN);				
		}

			}
	
		}
	}


int circ_coord(int x, int n){
	return (x>n)*(x%n) + (x<0)*((-x)%n);

}




void DrawSegments2d(segment *segs,compartment *cmps, int n,double min_v,double max_v,int cur_i,int w,int h,double scale,Vector2 v){
	
	Vector2 soma_pos = {segs->x,segs->y};
	for(int i=0;i<n;++i){
		segment *seg = segs + i ;
		compartment *cmp = cmps + i; 
		double clipd_v = cmp->v*(cmp->v < max_v && cmp->v > min_v) + max_v*(cmp->v > max_v) + min_v*(cmp->v < min_v);
		double ratio = (clipd_v-min_v)/(max_v-min_v);
		Vector2 seg_pos = (Vector2){seg->x + w/2 - soma_pos.x ,seg->y+h/2 - soma_pos.y};
		Color c = (Color){.r = (int)(255*ratio), .g = 0,.b = (int)255*(1-ratio),.a = 255};
	
		if(seg->id == 0){
			Vector2 soma_pos = (Vector2){seg_pos.x ,seg_pos.y};
			if(i == cur_i){
			c = BLACK;
		}


		} {
			segment *paren = seg->p;
			Vector2 paren_pos =  (Vector2){paren->x-soma_pos.x + w/2 ,paren->y-soma_pos.y+h/2};
			DrawLineV(seg_pos, paren_pos, c);
			if(i == cur_i){
			c = BLACK;
		}	}
	
		}
	}





char **filter_by_first_char(char **lines,char chr){
	int num_lines = 1;
	char **nu_lines = malloc(sizeof(char*));
	while(strcmp(*lines, "\0") != 0){
		if(**lines == chr){
			++lines;
		}else {
			nu_lines[num_lines-1] = *lines;
			++lines;
			++num_lines;
			nu_lines = realloc(nu_lines, num_lines*sizeof(char *));
			
		}


	}
	nu_lines[num_lines-1] = "\0";
	return nu_lines;
}

char **drop_k(char **lines,int k){
	int num_lines = 1;
	char **nu_lines = malloc(sizeof(char*));
	while(strcmp(*lines, "\0") != 0){
			nu_lines[num_lines-1] = *(lines)+k;
			++lines;
			++num_lines;
			nu_lines = realloc(nu_lines, num_lines*sizeof(char *));
			
	}
	nu_lines[num_lines-1] = "\0";
	return nu_lines;
}

ptr_w_len char_to_double_array(char **lines,char col_delim){	
	int N = 1;
	double *dub_arry = malloc(sizeof(double)*N);
	while(strcmp(*lines, "\0") != 0){	
		char **str_doubles = break_at_delim(*lines, col_delim);
		while(strcmp(*str_doubles, "\0") != 0){
			dub_arry[N-1] = atof(*str_doubles);
			++str_doubles;
			++N;
			dub_arry = realloc(dub_arry, sizeof(double)*N);
		}
	}
	return (ptr_w_len){.ptr = dub_arry,.len = N};

}

//
//
//
//
//ENTERING
//MODELLAND


// Simulation constants
#define DT 0.01    // Time step (ms)
#define T_MAX 100.0 // Total simulation time (ms)

// Physical constants
#define C_M 1.0   // Membrane capacitance (uF/cm^2)
#define G_L 0.3    // Leak conductance (mS/cm^2)
#define E_L -54.4  // Leak reversal potential (mV)

// Inline functions for gating calculations
static inline double alpha(double V, gating_params params) {
    return params.scale * (V - params.V_m) / (1 - exp(-(V - params.V_m) / params.k));
}

static inline double beta(double V, gating_params params) {
    return params.scale * exp(-(V - params.V_m) / params.k);
}

static inline double inf(double alpha, double beta) {
    return alpha / (alpha + beta);
}


// Simulation constants
// #define DT 0.01    // Time step (ms)
// #define T_MAX 100.0 // Total simulation time (ms)
//
// // Physical constants
// #define C_M 1.0   // Membrane capacitance (uF/cm^2)
// #define G_L 0.3    // Leak conductance (mS/cm^2)
// #define E_L -54.4  // Leak reversal potential (mV)
//
// // Struct for gating variable parameters
// typedef struct {
//     double V_m;    // Midpoint voltage (mV)
//         double k;      // Slope factor (mV)
//             double scale;  // Scaling factor
//             } gating_params;
//
//             // Inline functions for gating calculations
//             static inline double alpha(double V, gating_params params) {
//                 return params.scale * (V - params.V_m) / (1 - exp(-(V - params.V_m) / params.k));
//                 }
//
//                 static inline double beta(double V, gating_params params) {
//                     return params.scale * exp(-(V - params.V_m) / params.k);
//                     }
//
//                     static inline double inf(double alpha, double beta) {
//                         return alpha / (alpha + beta);
//                         }
//
//                         static inline double tau(double alpha, double beta) {
//                             return 1.0 / (alpha + beta);
//                             }
//

static inline double tau(double alpha, double beta) {
    return 1.0 / (alpha + beta);
}

// Struct for ion channel
typedef struct {
    int n_m;              // Number of gating variables
    double *m;            // Array of gating variable values
    double *tau;          // Array of time constants
    int *exponents;       // Array of exponents for each gating variable
    gating_params *m_params; // Array of gating variable parameters
    double (*m_inf)(double V, gating_params params); // Pointer to m_inf function
    double g_max;         // Maximum conductance (mS/cm^2)
    double reversal;      // Reversal potential (mV)
} Channel;

// Function to compute steady-state value
double compute_m_inf(double V, gating_params params) {
    double a = alpha(V, params);
    double b = beta(V, params);
    return inf(a, b);
}

// Function to initialize a channel
void initialize_channel(Channel *channel, int n_m, gating_params *params, int *exponents, double g_max, double reversal) {
    channel->n_m = n_m;
    channel->m = calloc(n_m, sizeof(double));
    channel->tau = calloc(n_m, sizeof(double));
    channel->exponents = malloc(n_m * sizeof(int));
    channel->m_params = params;
    channel->m_inf = compute_m_inf;
    channel->g_max = g_max;
    channel->reversal = reversal;

    // Initialize gating variables and exponents
    for (int i = 0; i < n_m; i++) {
        channel->m[i] = 0.0; // Initial gating variable values
        channel->tau[i] = 0.0;
        channel->exponents[i] = exponents[i];
    }
}

// Function to compute the channel current (without using `pow`)
double compute_channel_current(Channel *channel, double V) {
    double activation = 1.0;
    for (int i = 0; i < channel->n_m; i++) {
        // Compute steady-state and time constant
        double a = alpha(V, channel->m_params[i]);
        double b = beta(V, channel->m_params[i]);
        channel->tau[i] = tau(a, b);
        // Manual exponentiation instead of using `pow`
        double m_exponent = 1.0;
        for (int j = 0; j < channel->exponents[i]; j++) {
            m_exponent *= channel->m[i];
        }
        activation *= m_exponent;
        // Update gating variable
        channel->m[i] += DT * (channel->m_inf(V, channel->m_params[i]) - channel->m[i]) / channel->tau[i];
    }
    return channel->g_max * activation * (V - channel->reversal);
}


void upd_cmp(compartment *cmp){
	cmp->v = cmp->v + cmp->temp_v;
	cmp->temp_v = 0;
}

compartment *segments_to_compartnents(segment *segs,int n){
	compartment *cmp = malloc(n*sizeof(compartment));
	for(int i=0;i<n;++i){
		segment *c_seg = segs+i;
		(cmp+i)->d  = 0;
		while(c_seg->p_id != c_seg->id){
			(cmp+i)->d  += c_seg->l;
			c_seg = segs + c_seg->p_id;
		}
		//printf("i = %d",i);
		(cmp+i)->id = (segs+i)->id;
		(cmp+i)->n_ne = (segs+i)->n_c + 1;
		//printf(":: n_ne = %d \n",(cmp+i)->n_ne);
		(cmp+i)->id_ne = realloc((cmp+i)->id_ne,(cmp+i)->n_ne*sizeof(int));
		(cmp+i)->g_ne = realloc((cmp+i)->g_ne,(cmp+i)->n_ne*sizeof(double));
		(cmp+i)->id_ne[0] = (segs+i)->p_id;
		(cmp+i)->g_ne[0] = (segs+i)->g_p;
		(cmp+i)->S = (segs+i)->S;
		(cmp+i)->v = (segs+i)->v;
		(cmp+i)->temp_v = 0;
		(cmp+i)->type = (segs+i)->type;
		for(int i_ne = 1;i_ne<(cmp+i)->n_ne;++i_ne){
			(cmp+i)->g_ne[i_ne] = (segs+i)->g_c[i_ne-1];
			(cmp+i)->id_ne[i_ne] = (segs+i)->c_id[i_ne-1];
			
		}
	}
	return cmp;
}

void upd_net(net *net){
//#pragma omp parallel for
	for(int i =0;i<net->n;++i){

		compartment *cur = net->cmp + i;
		net->I_m(net->cmp + i,net->dt);
		net->I_l(cur,net->dt);
	}
//#pragma omp parallel for
	for(int i =0;i<net->n;++i){
		compartment *cur = net->cmp + i;

		net->upd_cmp(cur);	
	
	}
//#pragma omp parallel for
	for(int i =0;i<net->n;++i){
		compartment *cur = net->cmp + i;
		for(int i_chan = 0;i_chan < cur->n_chans;++i_chan){
		net->upd_chan(cur,cur->chans + i_chan,net->dt);	}
	}
}



double g_from_chan(compartment *cmp,channel *chan,double g_max){
	double g = g_max;
	for(int i = 0;i<chan->n_m;++i){
		for(int j = 0;j<chan->exponents[i];++j){
			g = g*chan->m[i];
		}	
	}
	return g;
}

void update_chan(compartment *cmp, channel *chan,double dt){
	for(int i_m = 0;i_m<chan->n_m;++i_m){
		chan->m[i_m] = (dt/chan->tau[i_m])*((chan->m_inf + i_m)(cmp->v) -chan->m[i_m]); 
		}
	}




void I_m(compartment *cmp,double dt){
	for(int i = 0;i<cmp->n_chans;++i){
		cmp->temp_v = cmp->temp_v + dt*g_from_chan(cmp,cmp->chans + i,cmp->g_max[i])*((cmp->chans + i)->reversal-cmp->v); 
	}
	
}

void I_l(compartment *cmp,double dt){
		compartment *cmp_0 = cmp - cmp->id;
		for(int j = 0;j<cmp->n_ne;++j){
			compartment *ne = cmp_0 + cmp->id_ne[j];
			cmp->temp_v = cmp->temp_v + dt*cmp->g_ne[j]*(ne->v - cmp->v);
		}
}


net *net_skeleton(net *net){
	net->I_l = &I_l;
	net->I_m = &I_m;
	net->upd_chan = &update_chan;
	net->upd = &upd_net;
	net->upd_cmp = &upd_cmp;
	net->n = 0;
	//net->dt = 0;
	return net;
}

void add_cmps(net *net,compartment *cmp,int n){
	net->n = n;
	net->cmp = cmp;
}

void add_chan(net *net,channel *chan,double (*g_max)(compartment *cmp),int (*cmp_filter)(compartment *cmp)){
	printf("we here chan net has %d",net->n);
	for(int i_cmp = 0; i_cmp < net->n;++i_cmp){
		printf("we here chan %d",i_cmp);
		compartment *cmp = net->cmp + i_cmp;
		if(1){
			cmp->chans = realloc(cmp->chans, (cmp->n_chans + 1)*sizeof(channel));
			cmp->chans[cmp->n_chans]  = *chan; 
			
			cmp->g_max = realloc(cmp->g_max, (cmp->n_chans + 1)*sizeof(channel));
			cmp->g_max[cmp->n_chans] = g_max(cmp);

			cmp->n_chans += 1;

		
		}
	}

}



double g_l(compartment *cmp){
	return G_L;
}

channel leak = {.exponents = 0,.n_m = 0,.reversal = 0 };


int all(compartment *cmp){
	return 1;
}


int main(){




FILE *ptr;

ptr = fopen("Acker2008.swc","r");

fseek(ptr,0L,SEEK_END);
long sz = ftell(ptr);
rewind(ptr);

char *str = malloc(sz+1);
fread(str, 1, sz, ptr);
str[sz] = '\0';

char **lines = drop_k(filter_by_first_char(break_at_delim(str, '\n'),'#'),1);
int n = 0;
while(strcmp(*(lines+n),"\0") != 0){
	++n;
}

point* segmns = malloc(n*sizeof(point));
init_segs(lines, n,segmns,0.7);

segment *segs = new_segments(segmns,n);
add_gs(segs,n, 0.01);
compartment *cmp = segments_to_compartnents(segs, n);



net *nt = malloc(sizeof(net));
nt = net_skeleton(nt);
add_cmps(nt, cmp, n);



double gl = 0.01;
printf("we hete \n");
add_chan(nt, &leak, &g_l,&all);
nt->dt = 0.001;

print_net(nt, 10);



printf("avail threads = %d \n",omp_get_max_threads());




// Initialization
//--------------------------------------------------------------------------------------
const int screenWidth = 800;
const int screenHeight = 450;

InitWindow(screenWidth, screenHeight,"");

int cur_i = 0;
int draw_every = 100, k = 0;
 // Main game loop
    while (!WindowShouldClose())        // Detect window close button or ESC key
    {
	cur_i = (cur_i > 0 && cur_i < n)*cur_i  + (cur_i > n)*(n-cur_i) + (cur_i <0)*(n+cur_i);
	
        //----------------------------------------------------------------------------------
	//
	if(k % draw_every == 0){
	if(IsKeyPressed(KEY_TAB)){
		cur_i = cur_i + 1;
	}
        BeginDrawing();

            ClearBackground(RAYWHITE);
	    	DrawSegments2d(segs, cmp,n, 100,  0, cur_i,screenWidth,screenHeight);
           	DrawText(TextFormat("fps = %d",GetFPS()), 0, 0, 20, BLACK);        
        	DrawText(TextFormat("cur_i = %d",cur_i), 0, 30, 20, BLACK);        
		EndDrawing();}
	++k;
	
	if(IsKeyDown(KEY_I)){
		(nt->cmp + cur_i)->v += nt->dt; 

	}	
	if(IsKeyDown(KEY_O)){
		(nt->cmp + cur_i)->v -= nt->dt; 

	}
	nt->upd(nt);
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

return 0;}

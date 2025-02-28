/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__capr
#define _nrn_initial _nrn_initial__capr
#define nrn_cur _nrn_cur__capr
#define _nrn_current _nrn_current__capr
#define nrn_jacob _nrn_jacob__capr
#define nrn_state _nrn_state__capr
#define _net_receive _net_receive__capr 
#define pmp pmp__capr 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define pump _p[0]
#define pump_columnindex 0
#define pumpca _p[1]
#define pumpca_columnindex 1
#define Dpump _p[2]
#define Dpump_columnindex 2
#define Dpumpca _p[3]
#define Dpumpca_columnindex 3
#define cai _p[4]
#define cai_columnindex 4
#define Dcai _p[5]
#define Dcai_columnindex 5
#define cao _p[6]
#define cao_columnindex 6
#define ica _p[7]
#define ica_columnindex 7
#define ipump _p[8]
#define ipump_columnindex 8
#define ipump_last _p[9]
#define ipump_last_columnindex 9
#define voli _p[10]
#define voli_columnindex 10
#define area1 _p[11]
#define area1_columnindex 11
#define c1 _p[12]
#define c1_columnindex 12
#define c2 _p[13]
#define c2_columnindex 13
#define c3 _p[14]
#define c3_columnindex 14
#define c4 _p[15]
#define c4_columnindex 15
#define v _p[16]
#define v_columnindex 16
#define _g _p[17]
#define _g_columnindex 17
#define _ion_cao	*_ppvar[0]._pval
#define _ion_cai	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
#define _style_ca	*((int*)_ppvar[4]._pvoid)
#define diam	*_ppvar[5]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_capr", _hoc_setdata,
 0, 0
};
 /* declare global and static user variables */
#define car car_capr
 double car = 5e-05;
#define k4 k4_capr
 double k4 = 5;
#define k3 k3_capr
 double k3 = 500;
#define k2 k2_capr
 double k2 = 250000;
#define k1 k1_capr
 double k1 = 5e+08;
#define pumpdens pumpdens_capr
 double pumpdens = 3e-14;
#define tau tau_capr
 double tau = 1e+09;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "car_capr", "mM",
 "tau_capr", "ms",
 "k1_capr", "/mM-s",
 "k2_capr", "/s",
 "k3_capr", "/s",
 "k4_capr", "/mM-s",
 "pumpdens_capr", "mol/cm2",
 "pump_capr", "mol/cm2",
 "pumpca_capr", "mol/cm2",
 0,0
};
 static double cai0 = 0;
 static double delta_t = 0.01;
 static double pumpca0 = 0;
 static double pump0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "car_capr", &car_capr,
 "tau_capr", &tau_capr,
 "k1_capr", &k1_capr,
 "k2_capr", &k2_capr,
 "k3_capr", &k3_capr,
 "k4_capr", &k4_capr,
 "pumpdens_capr", &pumpdens_capr,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"capr",
 0,
 0,
 "pump_capr",
 "pumpca_capr",
 0,
 0};
 static Symbol* _morphology_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 18, _prop);
 	/*initialize range parameters*/
 	_prop->param = _p;
 	_prop->param_size = 18;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_morphology_sym);
 	_ppvar[5]._pval = &prop_ion->param[0]; /* diam */
 prop_ion = need_memb(_ca_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[0]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[1]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 	_ppvar[4]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for ca */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _capmp_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	_morphology_sym = hoc_lookup("morphology");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 3);
  _extcall_thread = (Datum*)ecalloc(2, sizeof(Datum));
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 18, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "#ca_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
  hoc_register_dparam_semantics(_mechtype, 5, "diam");
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 capr /home/marko/magy/thesis/sim/capmp.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define PI _nrnunit_PI[_nrnunit_use_legacy_]
static double _nrnunit_PI[2] = {0x1.921fb54442d18p+1, 3.14159}; /* 3.14159265358979312 */
 
#define FARADAY _nrnunit_FARADAY[_nrnunit_use_legacy_]
static double _nrnunit_FARADAY[2] = {0x1.78e555060882cp+16, 96485.3}; /* 96485.3321233100141 */
 static double volo = 1;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 extern double *_nrn_thread_getelm(SparseObj*, int, int);
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  0
 static int _spth1 = 1;
 static int _cvspth1 = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3]; static double *_temp1;
 static int pmp();
 
static int pmp (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<3;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
}  
_RHS1(1) *= ( voli) ;
_MATELM1(1, 1) *= ( voli); 
_RHS1(2) *= ( ( 1e10 ) * area1) ;
_MATELM1(2, 2) *= ( ( 1e10 ) * area1);  }
 /* COMPARTMENT voli {
     cai }
   */
 /* COMPARTMENT ( 1e10 ) * area1 {
     pump pumpca }
   */
 /* COMPARTMENT volo * ( 1e15 ) {
     }
   */
 /* ~ car <-> cai ( 1.0 / tau , 1.0 / tau )*/
 f_flux =  1.0 / tau * car ;
 b_flux =  1.0 / tau * cai ;
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  1.0 / tau ;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ cai + pump <-> pumpca ( c1 , c2 )*/
 f_flux =  c1 * pump * cai ;
 b_flux =  c2 * pumpca ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 1) -= (f_flux - b_flux);
 
 _term =  c1 * cai ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  += _term;
 _term =  c1 * pump ;
 _MATELM1( 2 ,1)  += _term;
 _MATELM1( 1 ,1)  += _term;
 _term =  c2 ;
 _MATELM1( 2 ,0)  -= _term;
 _MATELM1( 1 ,0)  -= _term;
 /*REACTION*/
  /* ~ pumpca <-> pump + cao ( c3 , c4 )*/
 f_flux =  c3 * pumpca ;
 b_flux =  c4 * cao * pump ;
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  c3 ;
 _MATELM1( 2 ,0)  -= _term;
 _term =  c4 * cao ;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  ipump = ( 1e-4 ) * 2.0 * FARADAY * ( f_flux - b_flux ) / area1 ;
   /* ~ cai < < ( - ( ica - ipump_last ) * area1 / ( 2.0 * FARADAY ) * ( 1e4 ) )*/
 f_flux = b_flux = 0.;
 _RHS1( 1) += (b_flux =   ( - ( ica - ipump_last ) * area1 / ( 2.0 * FARADAY ) * ( 1e4 ) ) );
 /*FLUX*/
   /* pump + pumpca = ( 1e10 ) * area1 * pumpdens */
 _RHS1(0) =  ( 1e10 ) * area1 * pumpdens;
 _MATELM1(0, 0) = 1 * ( ( 1e10 ) * area1);
 _RHS1(0) -= pumpca * ( ( 1e10 ) * area1) ;
 _MATELM1(0, 2) = 1 * ( ( 1e10 ) * area1);
 _RHS1(0) -= pump * ( ( 1e10 ) * area1) ;
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<3;_i++) _p[_dlist1[_i]] = 0.0;}
 /* COMPARTMENT voli {
   cai }
 */
 /* COMPARTMENT ( 1e10 ) * area1 {
   pump pumpca }
 */
 /* COMPARTMENT volo * ( 1e15 ) {
   }
 */
 /* ~ car <-> cai ( 1.0 / tau , 1.0 / tau )*/
 f_flux =  1.0 / tau * car ;
 b_flux =  1.0 / tau * cai ;
 Dcai += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ cai + pump <-> pumpca ( c1 , c2 )*/
 f_flux =  c1 * pump * cai ;
 b_flux =  c2 * pumpca ;
 Dpump -= (f_flux - b_flux);
 Dcai -= (f_flux - b_flux);
 Dpumpca += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ pumpca <-> pump + cao ( c3 , c4 )*/
 f_flux =  c3 * pumpca ;
 b_flux =  c4 * cao * pump ;
 Dpumpca -= (f_flux - b_flux);
 Dpump += (f_flux - b_flux);
 
 /*REACTION*/
  ipump = ( 1e-4 ) * 2.0 * FARADAY * ( f_flux - b_flux ) / area1 ;
 /* ~ cai < < ( - ( ica - ipump_last ) * area1 / ( 2.0 * FARADAY ) * ( 1e4 ) )*/
 f_flux = b_flux = 0.;
 Dcai += (b_flux =   ( - ( ica - ipump_last ) * area1 / ( 2.0 * FARADAY ) * ( 1e4 ) ) );
 /*FLUX*/
   /* pump + pumpca = ( 1e10 ) * area1 * pumpdens */
 /*CONSERVATION*/
 _p[_dlist1[0]] /= ( ( 1e10 ) * area1);
 _p[_dlist1[1]] /= ( voli);
 _p[_dlist1[2]] /= ( ( 1e10 ) * area1);
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<3;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
}  
_RHS1(0) *= ( ( 1e10 ) * area1) ;
_MATELM1(0, 0) *= ( ( 1e10 ) * area1); 
_RHS1(1) *= ( voli) ;
_MATELM1(1, 1) *= ( voli); 
_RHS1(2) *= ( ( 1e10 ) * area1) ;
_MATELM1(2, 2) *= ( ( 1e10 ) * area1);  }
 /* COMPARTMENT voli {
 cai }
 */
 /* COMPARTMENT ( 1e10 ) * area1 {
 pump pumpca }
 */
 /* COMPARTMENT volo * ( 1e15 ) {
 }
 */
 /* ~ car <-> cai ( 1.0 / tau , 1.0 / tau )*/
 _term =  1.0 / tau ;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ cai + pump <-> pumpca ( c1 , c2 )*/
 _term =  c1 * cai ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  += _term;
 _MATELM1( 0 ,2)  -= _term;
 _term =  c1 * pump ;
 _MATELM1( 2 ,1)  += _term;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  c2 ;
 _MATELM1( 2 ,0)  -= _term;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ pumpca <-> pump + cao ( c3 , c4 )*/
 _term =  c3 ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 2 ,0)  -= _term;
 _term =  c4 * cao ;
 _MATELM1( 0 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /* ~ cai < < ( - ( ica - ipump_last ) * area1 / ( 2.0 * FARADAY ) * ( 1e4 ) )*/
 /*FLUX*/
   /* pump + pumpca = ( 1e10 ) * area1 * pumpdens */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cao = _ion_cao;
  cai = _ion_cai;
  cai = _ion_cai;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  _ion_cai = cai;
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 	_pv[1] = &(_ion_cai);
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 3, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cao = _ion_cao;
  cai = _ion_cai;
  cai = _ion_cai;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(_thread[_cvspth1]._pvoid);
   _nrn_destroy_sparseobj_thread(_thread[_spth1]._pvoid);
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  pumpca = pumpca0;
  pump = pump0;
 {
   voli = PI * pow( ( diam / 2.0 ) , 2.0 ) * 1.0 ;
   area1 = 2.0 * PI * ( diam / 2.0 ) * 1.0 ;
   c1 = ( 1e7 ) * area1 * k1 ;
   c2 = ( 1e7 ) * area1 * k2 ;
   c3 = ( 1e7 ) * area1 * k3 ;
   c4 = ( 1e7 ) * area1 * k4 ;
    _ss_sparse_thread(&_thread[_spth1]._pvoid, 3, _slist1, _dlist1, _p, &t, dt, pmp, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cao = _ion_cao;
  cai = _ion_cai;
  cai = _ion_cai;
 initmodel(_p, _ppvar, _thread, _nt);
  _ion_cai = cai;
   nrn_wrote_conc(_ca_sym, (&(_ion_cai)) - 1, _style_ca);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ipump_last = ipump ;
   ica = ipump ;
   }
 _current += ica;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  cao = _ion_cao;
  cai = _ion_cai;
  cai = _ion_cai;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_cai = cai;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  cao = _ion_cao;
  cai = _ion_cai;
  cai = _ion_cai;
 {  sparse_thread(&_thread[_spth1]._pvoid, 3, _slist1, _dlist1, _p, &t, dt, pmp, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 3; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 }  _ion_cai = cai;
 }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = pumpca_columnindex;  _dlist1[0] = Dpumpca_columnindex;
 _slist1[1] = cai_columnindex;  _dlist1[1] = Dcai_columnindex;
 _slist1[2] = pump_columnindex;  _dlist1[2] = Dpump_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/marko/magy/thesis/sim/capmp.mod";
static const char* nmodl_file_text = 
  ":  capump.mod plus a \"reservoir\" used to initialize cai to desired concentrations\n"
  "\n"
  "UNITS {\n"
  "	(mM) = (milli/liter)\n"
  "	(mA) = (milliamp)\n"
  "	(um) = (micron)\n"
  "	(mol) = (1)\n"
  "	PI = (pi) (1)\n"
  "	FARADAY = (faraday) (coulomb)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX capr\n"
  "	USEION ca READ cao, cai WRITE cai, ica\n"
  "	GLOBAL k1, k2, k3, k4\n"
  "	GLOBAL car, tau\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	pump	(mol/cm2)\n"
  "	pumpca	(mol/cm2)\n"
  "	cai	(mM)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	car = 5e-5 (mM) : ca in reservoir, used to initialize cai to desired concentrations\n"
  "	tau = 1e9 (ms) : rate of equilibration between cai and car\n"
  "\n"
  "	k1 = 5e8	(/mM-s)\n"
  "	k2 = .25e6	(/s)\n"
  "	k3 = .5e3	(/s)\n"
  "	k4 = 5e0	(/mM-s)\n"
  "\n"
  "	pumpdens = 3e-14 (mol/cm2)\n"
  "}\n"
  "\n"
  "CONSTANT {\n"
  "	volo = 1 (liter)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	diam	(um)\n"
  "	cao	(mM)\n"
  "\n"
  "	ica (mA/cm2)\n"
  "	ipump (mA/cm2)\n"
  "	ipump_last (mA/cm2)\n"
  "	voli	(um3)\n"
  "	area1	(um2)\n"
  "	c1	(1+8 um5/ms)\n"
  "	c2	(1-10 um2/ms)\n"
  "	c3	(1-10 um2/ms)\n"
  "	c4	(1+8 um5/ms)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE pmp METHOD sparse\n"
  "	ipump_last = ipump\n"
  "	ica = ipump\n"
  "}\n"
  "\n"
  "KINETIC pmp {\n"
  "	COMPARTMENT voli {cai}\n"
  "	COMPARTMENT (1e10)*area1 {pump pumpca}\n"
  "	COMPARTMENT volo*(1e15) {cao car}\n"
  "\n"
  "	~ car <-> cai		(1(um3)/tau,1(um3)/tau)\n"
  "	~ cai + pump <-> pumpca		(c1,c2)\n"
  "	~ pumpca     <-> pump + cao	(c3,c4)\n"
  "\n"
  "	: note that forward flux here is the outward flux\n"
  "	ipump = (1e-4)*2*FARADAY*(f_flux - b_flux)/area1\n"
  "\n"
  "        : ipump_last vs ipump needed because of STEADYSTATE calculation\n"
  "        ~ cai << (-(ica - ipump_last)*area1/(2*FARADAY)*(1e4))\n"
  "\n"
  "	CONSERVE pump + pumpca = (1e10)*area1*pumpdens\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	:cylindrical coordinates; actually vol and area1/unit length\n"
  "	voli = PI*(diam/2)^2 * 1(um)\n"
  "	area1 = 2*PI*(diam/2) * 1(um)\n"
  "	c1 = (1e7)*area1 * k1\n"
  "	c2 = (1e7)*area1 * k2\n"
  "	c3 = (1e7)*area1 * k3\n"
  "	c4 = (1e7)*area1 * k4\n"
  "\n"
  "	SOLVE pmp STEADYSTATE sparse\n"
  "}\n"
  ;
#endif

/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib.h"
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
 
#define nrn_init _nrn_init__gaines_internode_flutstin
#define _nrn_initial _nrn_initial__gaines_internode_flutstin
#define nrn_cur _nrn_cur__gaines_internode_flutstin
#define _nrn_current _nrn_current__gaines_internode_flutstin
#define nrn_jacob _nrn_jacob__gaines_internode_flutstin
#define nrn_state _nrn_state__gaines_internode_flutstin
#define _net_receive _net_receive__gaines_internode_flutstin 
#define evaluate_fct evaluate_fct__gaines_internode_flutstin 
#define states states__gaines_internode_flutstin 
 
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
#define gkbar _p[0]
#define gkfbar _p[1]
#define gl _p[2]
#define ghcnbar _p[3]
#define ek _p[4]
#define el _p[5]
#define eq _p[6]
#define ik _p[7]
#define ikf _p[8]
#define ihcn _p[9]
#define il _p[10]
#define n_inf _p[11]
#define s_inf _p[12]
#define q_inf _p[13]
#define tau_n _p[14]
#define tau_s _p[15]
#define tau_q _p[16]
#define s _p[17]
#define n _p[18]
#define q _p[19]
#define Ds _p[20]
#define Dn _p[21]
#define Dq _p[22]
#define q10 _p[23]
#define v _p[24]
#define _g _p[25]
 
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
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_Exp(void);
 static void _hoc_evaluate_fct(void);
 static void _hoc_vtrap6(void);
 static void _hoc_vtrap5(void);
 static void _hoc_vtrap4(void);
 static void _hoc_vtrap3(void);
 static void _hoc_vtrap2(void);
 static void _hoc_vtrap1(void);
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
 "setdata_gaines_internode_flutstin", _hoc_setdata,
 "Exp_gaines_internode_flutstin", _hoc_Exp,
 "evaluate_fct_gaines_internode_flutstin", _hoc_evaluate_fct,
 "vtrap6_gaines_internode_flutstin", _hoc_vtrap6,
 "vtrap5_gaines_internode_flutstin", _hoc_vtrap5,
 "vtrap4_gaines_internode_flutstin", _hoc_vtrap4,
 "vtrap3_gaines_internode_flutstin", _hoc_vtrap3,
 "vtrap2_gaines_internode_flutstin", _hoc_vtrap2,
 "vtrap1_gaines_internode_flutstin", _hoc_vtrap1,
 0, 0
};
#define Exp Exp_gaines_internode_flutstin
#define vtrap6 vtrap6_gaines_internode_flutstin
#define vtrap5 vtrap5_gaines_internode_flutstin
#define vtrap4 vtrap4_gaines_internode_flutstin
#define vtrap3 vtrap3_gaines_internode_flutstin
#define vtrap2 vtrap2_gaines_internode_flutstin
#define vtrap1 vtrap1_gaines_internode_flutstin
 extern double Exp( _threadargsprotocomma_ double );
 extern double vtrap6( _threadargsprotocomma_ double );
 extern double vtrap5( _threadargsprotocomma_ double );
 extern double vtrap4( _threadargsprotocomma_ double );
 extern double vtrap3( _threadargsprotocomma_ double );
 extern double vtrap2( _threadargsprotocomma_ double );
 extern double vtrap1( _threadargsprotocomma_ double );
 /* declare global and static user variables */
#define aqC aqC_gaines_internode_flutstin
 double aqC = -12.2;
#define aqB aqB_gaines_internode_flutstin
 double aqB = -107.3;
#define aqA aqA_gaines_internode_flutstin
 double aqA = 0.00522;
#define asC asC_gaines_internode_flutstin
 double asC = -5;
#define asB asB_gaines_internode_flutstin
 double asB = -27;
#define asA asA_gaines_internode_flutstin
 double asA = 0.3;
#define anC anC_gaines_internode_flutstin
 double anC = 1.1;
#define anB anB_gaines_internode_flutstin
 double anB = -83.2;
#define anA anA_gaines_internode_flutstin
 double anA = 0.0462;
#define bqC bqC_gaines_internode_flutstin
 double bqC = -12.2;
#define bqB bqB_gaines_internode_flutstin
 double bqB = -107.3;
#define bqA bqA_gaines_internode_flutstin
 double bqA = 0.00522;
#define bsC bsC_gaines_internode_flutstin
 double bsC = -1;
#define bsB bsB_gaines_internode_flutstin
 double bsB = 10;
#define bsA bsA_gaines_internode_flutstin
 double bsA = 0.03;
#define bnC bnC_gaines_internode_flutstin
 double bnC = 10.5;
#define bnB bnB_gaines_internode_flutstin
 double bnB = -66;
#define bnA bnA_gaines_internode_flutstin
 double bnA = 0.0824;
#define vtraub vtraub_gaines_internode_flutstin
 double vtraub = -80;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gkbar_gaines_internode_flutstin", "mho/cm2",
 "gkfbar_gaines_internode_flutstin", "mho/cm2",
 "gl_gaines_internode_flutstin", "mho/cm2",
 "ek_gaines_internode_flutstin", "mV",
 "el_gaines_internode_flutstin", "mV",
 "eq_gaines_internode_flutstin", "mV",
 "ik_gaines_internode_flutstin", "mA/cm2",
 "ikf_gaines_internode_flutstin", "mA/cm2",
 "ihcn_gaines_internode_flutstin", "mA/cm2",
 "il_gaines_internode_flutstin", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double n0 = 0;
 static double q0 = 0;
 static double s0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vtraub_gaines_internode_flutstin", &vtraub_gaines_internode_flutstin,
 "anA_gaines_internode_flutstin", &anA_gaines_internode_flutstin,
 "anB_gaines_internode_flutstin", &anB_gaines_internode_flutstin,
 "anC_gaines_internode_flutstin", &anC_gaines_internode_flutstin,
 "bnA_gaines_internode_flutstin", &bnA_gaines_internode_flutstin,
 "bnB_gaines_internode_flutstin", &bnB_gaines_internode_flutstin,
 "bnC_gaines_internode_flutstin", &bnC_gaines_internode_flutstin,
 "asA_gaines_internode_flutstin", &asA_gaines_internode_flutstin,
 "asB_gaines_internode_flutstin", &asB_gaines_internode_flutstin,
 "asC_gaines_internode_flutstin", &asC_gaines_internode_flutstin,
 "bsA_gaines_internode_flutstin", &bsA_gaines_internode_flutstin,
 "bsB_gaines_internode_flutstin", &bsB_gaines_internode_flutstin,
 "bsC_gaines_internode_flutstin", &bsC_gaines_internode_flutstin,
 "aqA_gaines_internode_flutstin", &aqA_gaines_internode_flutstin,
 "aqB_gaines_internode_flutstin", &aqB_gaines_internode_flutstin,
 "aqC_gaines_internode_flutstin", &aqC_gaines_internode_flutstin,
 "bqA_gaines_internode_flutstin", &bqA_gaines_internode_flutstin,
 "bqB_gaines_internode_flutstin", &bqB_gaines_internode_flutstin,
 "bqC_gaines_internode_flutstin", &bqC_gaines_internode_flutstin,
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
 
#define _cvode_ieq _ppvar[0]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"gaines_internode_flutstin",
 "gkbar_gaines_internode_flutstin",
 "gkfbar_gaines_internode_flutstin",
 "gl_gaines_internode_flutstin",
 "ghcnbar_gaines_internode_flutstin",
 "ek_gaines_internode_flutstin",
 "el_gaines_internode_flutstin",
 "eq_gaines_internode_flutstin",
 0,
 "ik_gaines_internode_flutstin",
 "ikf_gaines_internode_flutstin",
 "ihcn_gaines_internode_flutstin",
 "il_gaines_internode_flutstin",
 "n_inf_gaines_internode_flutstin",
 "s_inf_gaines_internode_flutstin",
 "q_inf_gaines_internode_flutstin",
 "tau_n_gaines_internode_flutstin",
 "tau_s_gaines_internode_flutstin",
 "tau_q_gaines_internode_flutstin",
 0,
 "s_gaines_internode_flutstin",
 "n_gaines_internode_flutstin",
 "q_gaines_internode_flutstin",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 26, _prop);
 	/*initialize range parameters*/
 	gkbar = 0.002581;
 	gkfbar = 0.002568;
 	gl = 0.0002;
 	ghcnbar = 0.002232;
 	ek = -90;
 	el = -80;
 	eq = -54.9;
 	_prop->param = _p;
 	_prop->param_size = 26;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 1, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _gaines_internode_FLUTSTIN_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 26, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 gaines_internode_flutstin /Users/eliefarah/ascent/MOD_Files/gaines_internode_FLUTSTIN.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "Gaines Motor Axon Internode channels";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int evaluate_fct(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   evaluate_fct ( _threadargscomma_ v ) ;
   Ds = ( s_inf - s ) / tau_s ;
   Dn = ( n_inf - n ) / tau_n ;
   Dq = ( q_inf - q ) / tau_q ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 evaluate_fct ( _threadargscomma_ v ) ;
 Ds = Ds  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_s )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_n )) ;
 Dq = Dq  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_q )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   evaluate_fct ( _threadargscomma_ v ) ;
    s = s + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_s)))*(- ( ( ( s_inf ) ) / tau_s ) / ( ( ( ( - 1.0 ) ) ) / tau_s ) - s) ;
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_n)))*(- ( ( ( n_inf ) ) / tau_n ) / ( ( ( ( - 1.0 ) ) ) / tau_n ) - n) ;
    q = q + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_q)))*(- ( ( ( q_inf ) ) / tau_q ) / ( ( ( ( - 1.0 ) ) ) / tau_q ) - q) ;
   }
  return 0;
}
 
static int  evaluate_fct ( _threadargsprotocomma_ double _lv ) {
   double _la , _lb , _lv2 ;
 _la = q10 * vtrap1 ( _threadargscomma_ _lv ) ;
   _lb = q10 * vtrap2 ( _threadargscomma_ _lv ) ;
   tau_s = 1.0 / ( _la + _lb ) ;
   s_inf = _la / ( _la + _lb ) ;
   _la = q10 * vtrap3 ( _threadargscomma_ _lv ) ;
   _lb = q10 * vtrap4 ( _threadargscomma_ _lv ) ;
   tau_n = 1.0 / ( _la + _lb ) ;
   n_inf = _la / ( _la + _lb ) ;
   _la = q10 * vtrap5 ( _threadargscomma_ _lv ) ;
   _lb = q10 * vtrap6 ( _threadargscomma_ _lv ) ;
   tau_q = 1.0 / ( _la + _lb ) ;
   q_inf = _la / ( _la + _lb ) ;
   _lv2 = _lv - vtraub ;
    return 0; }
 
static void _hoc_evaluate_fct(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 evaluate_fct ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap1 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap1;
 _lvtrap1 = asA / ( Exp ( _threadargscomma_ ( _lx + asB ) / asC ) + 1.0 ) ;
   
return _lvtrap1;
 }
 
static void _hoc_vtrap1(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap1 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap2 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap2;
 _lvtrap2 = bsA / ( Exp ( _threadargscomma_ ( _lx + bsB ) / bsC ) + 1.0 ) ;
   
return _lvtrap2;
 }
 
static void _hoc_vtrap2(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap2 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap3 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap3;
 _lvtrap3 = anA * ( _lx - anB ) / ( 1.0 - Exp ( _threadargscomma_ ( anB - _lx ) / anC ) ) ;
   
return _lvtrap3;
 }
 
static void _hoc_vtrap3(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap3 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap4 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap4;
 _lvtrap4 = bnA * ( bnB - _lx ) / ( 1.0 - Exp ( _threadargscomma_ ( _lx - bnB ) / bnC ) ) ;
   
return _lvtrap4;
 }
 
static void _hoc_vtrap4(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap4 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap5 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap5;
 _lvtrap5 = aqA * Exp ( _threadargscomma_ ( _lx - aqB ) / aqC ) ;
   
return _lvtrap5;
 }
 
static void _hoc_vtrap5(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap5 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap6 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap6;
 _lvtrap6 = bqA / Exp ( _threadargscomma_ ( _lx - bqB ) / bqC ) ;
   
return _lvtrap6;
 }
 
static void _hoc_vtrap6(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap6 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double Exp ( _threadargsprotocomma_ double _lx ) {
   double _lExp;
 if ( _lx < - 100.0 ) {
     _lExp = 0.0 ;
     }
   else {
     _lExp = exp ( _lx ) ;
     }
   
return _lExp;
 }
 
static void _hoc_Exp(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  Exp ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
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
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
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
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  n = n0;
  q = q0;
  s = s0;
 {
   double _lt_howells ;
 _lt_howells = 34.0 ;
   q10 = pow( 3.0 , ( ( celsius - _lt_howells ) / 10.0 ) ) ;
   evaluate_fct ( _threadargscomma_ v ) ;
   s = s_inf ;
   n = n_inf ;
   q = q_inf ;
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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ik = gkbar * s * ( v - ek ) ;
   ikf = gkfbar * n * n * n * n * ( v - ek ) ;
   ihcn = ghcnbar * q * ( v - eq ) ;
   il = gl * ( v - el ) ;
   }
 _current += ikf;
 _current += ik;
 _current += ihcn;
 _current += il;

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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
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
 {   states(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(s) - _p;  _dlist1[0] = &(Ds) - _p;
 _slist1[1] = &(n) - _p;  _dlist1[1] = &(Dn) - _p;
 _slist1[2] = &(q) - _p;  _dlist1[2] = &(Dq) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/eliefarah/ascent/MOD_Files/gaines_internode_FLUTSTIN.mod";
static const char* nmodl_file_text = 
  "TITLE Gaines Motor Axon Internode channels\n"
  "\n"
  ": 2/02\n"
  ": Cameron C. McIntyre\n"
  ":\n"
  ": Fast Na+, Persistant Na+, Slow K+, and Leakage currents \n"
  ": responsible for nodal action potential\n"
  ": Iterative equations H-H notation rest = -80 mV\n"
  ":\n"
  ": This model is described in detail in:\n"
  ":\n"
  ": McIntyre CC, Richardson AG, and Grill WM. Modeling the excitability of\n"
  ": mammalian nerve fibers: influence of afterpotentials on the recovery\n"
  ": cycle. Journal of Neurophysiology 87:995-1006, 2002.\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX gaines_internode_flutstin\n"
  "	NONSPECIFIC_CURRENT ikf\n"
  "	NONSPECIFIC_CURRENT ik : Slow potassium\n"
  "	NONSPECIFIC_CURRENT ihcn : HCN channel\n"
  "	NONSPECIFIC_CURRENT il\n"
  "	RANGE gkbar, gl, gkfbar, ghcnbar, ena, ek, el, eq\n"
  "	RANGE s_inf, n_inf, q_inf	: gating variables, s -> slow potassium, n -> fast potassium, q -> HCN\n"
  "	RANGE tau_s, tau_n, tau_q\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "\n"
  "	gkbar   = 0.002581 	(mho/cm2)\n"
  "	gkfbar = 0.002568	(mho/cm2) : Different between flut and mysa/stin\n"
  "	gl	= 0.0002 (mho/cm2)	: Different between flut and mysa/stin\n"
  "	ghcnbar = 0.002232\n"
  "	ek  = -90.0 (mV)\n"
  "	el	= -80.0 (mV)\n"
  "	eq = -54.9	(mV)\n"
  "	celsius		(degC)\n"
  "	dt              (ms)\n"
  "	v               (mV)\n"
  "	vtraub=-80\n"
  "\n"
  "	anA = 0.0462\n"
  "	anB = -83.2\n"
  "	anC = 1.1\n"
  "	bnA = 0.0824\n"
  "	bnB = -66\n"
  "	bnC = 10.5\n"
  "	asA = 0.3\n"
  "	asB = -27\n"
  "	asC = -5\n"
  "	bsA = 0.03\n"
  "	bsB = 10\n"
  "	bsC = -1\n"
  "	aqA = 0.00522\n"
  "	aqB = -107.3\n"
  "	aqC = -12.2\n"
  "	bqA = 0.00522\n"
  "	bqB = -107.3\n"
  "	bqC = -12.2\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	s n q\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	ik      (mA/cm2)\n"
  "	ikf		(mA/cm2)\n"
  "	ihcn	(mA/cm2)\n"
  "	il      (mA/cm2)\n"
  "	n_inf\n"
  "	s_inf\n"
  "	q_inf\n"
  "	tau_n\n"
  "	tau_s\n"
  "	tau_q\n"
  "	q10\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "	ik   = gkbar * s * (v - ek)\n"
  "	ikf = gkfbar * n*n*n*n * (v-ek)\n"
  "	ihcn = ghcnbar * q * (v-eq)\n"
  "	il   = gl * (v - el)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {   : exact Hodgkin-Huxley equations\n"
  "    evaluate_fct(v)\n"
  "	s' = (s_inf - s) / tau_s\n"
  "	n' = (n_inf - n) / tau_n\n"
  "	q' = (q_inf - q) / tau_q\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL { LOCAL t_howells\n"
  ":\n"
  ":	Q10 adjustment\n"
  ":	According to Howells 2012, s n q have Q10 of 3\n"
  ":\n"
  "	t_howells = 34 :Temperature for Howells measurements, for Q10 conversion\n"
  "\n"
  "	q10 = 3.0 ^ ((celsius-t_howells)/ 10 )\n"
  "\n"
  "	evaluate_fct(v)\n"
  "	s = s_inf\n"
  "	n = n_inf\n"
  "	q = q_inf\n"
  "}\n"
  "\n"
  "PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2\n"
  "\n"
  "	a = q10*vtrap1(v)\n"
  "	b = q10*vtrap2(v)\n"
  "	tau_s = 1 / (a + b)\n"
  "	s_inf = a / (a + b)\n"
  "\n"
  "	a = q10*vtrap3(v)\n"
  "	b = q10*vtrap4(v)\n"
  "	tau_n = 1 / (a + b)\n"
  "	n_inf = a / (a + b)\n"
  "\n"
  "	a = q10*vtrap5(v)\n"
  "	b = q10*vtrap6(v)\n"
  "	tau_q = 1 / (a + b)\n"
  "	q_inf = a / (a + b)\n"
  "\n"
  "	v2 = v - vtraub : convert to traub convention\n"
  "}\n"
  "\n"
  ":FUNCTION vtrap(x) {\n"
  ":	if (x < -50) {\n"
  ":		vtrap = 0\n"
  ":	}else{\n"
  ":		vtrap = bsA / (Exp((x+bsB)/bsC) + 1)\n"
  ":	}\n"
  ":}\n"
  "\n"
  "FUNCTION vtrap1(x) {\n"
  "	:Alpha s gating\n"
  "	vtrap1 = asA/(Exp((x+asB)/asC) + 1)\n"
  "}\n"
  "\n"
  "FUNCTION vtrap2(x) {\n"
  "	:Beta s gating\n"
  "	vtrap2 = bsA/(Exp((x+bsB)/bsC) + 1)\n"
  "}\n"
  "\n"
  "FUNCTION vtrap3(x) {\n"
  "	:Alpha n gating\n"
  "	vtrap3 = anA * (x-anB) / (1 - Exp((anB-x)/anC))\n"
  "}\n"
  "\n"
  "FUNCTION vtrap4(x) { \n"
  "	:Beta n gating\n"
  "	vtrap4 = bnA * (bnB-x) / (1 - Exp((x-bnB)/bnC))\n"
  "}\n"
  "\n"
  "FUNCTION vtrap5(x) { \n"
  "	:Alpha hcn (q) gating\n"
  "	vtrap5 = aqA * Exp((x-aqB)/aqC)\n"
  "}\n"
  "\n"
  "FUNCTION vtrap6(x) {\n"
  "	:Beta hcn (q) gating\n"
  "	vtrap6 = bqA / Exp((x-bqB)/bqC)\n"
  "}\n"
  "\n"
  "FUNCTION Exp(x) {\n"
  "	if (x < -100) {\n"
  "		Exp = 0\n"
  "	}else{\n"
  "		Exp = exp(x)\n"
  "	}\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif

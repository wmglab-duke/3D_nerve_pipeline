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
 
#define nrn_init _nrn_init__SK
#define _nrn_initial _nrn_initial__SK
#define nrn_cur _nrn_cur__SK
#define _nrn_current _nrn_current__SK
#define nrn_jacob _nrn_jacob__SK
#define nrn_state _nrn_state__SK
#define _net_receive _net_receive__SK 
#define state state__SK 
 
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
#define gskbar _p[0]
#define sstau _p[1]
#define cal _p[2]
#define cascale _p[3]
#define gsk _p[4]
#define ca_i _p[5]
#define q _p[6]
#define ek _p[7]
#define ik _p[8]
#define ica _p[9]
#define Dca_i _p[10]
#define Dq _p[11]
#define v _p[12]
#define _g _p[13]
#define _ion_ica	*_ppvar[0]._pval
#define _ion_ek	*_ppvar[1]._pval
#define _ion_ik	*_ppvar[2]._pval
#define _ion_dikdv	*_ppvar[3]._pval
#define area	*_ppvar[4]._pval
 
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
 static void _hoc_alphaq(void);
 static void _hoc_betaq(void);
 static void _hoc_exp1(void);
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
 "setdata_SK", _hoc_setdata,
 "alphaq_SK", _hoc_alphaq,
 "betaq_SK", _hoc_betaq,
 "exp1_SK", _hoc_exp1,
 0, 0
};
#define alphaq alphaq_SK
#define betaq betaq_SK
#define exp1 exp1_SK
 extern double alphaq( _threadargsprotocomma_ double );
 extern double betaq( _threadargsprotocomma_ double );
 extern double exp1( _threadargsprotocomma_ double , double , double , double );
 /* declare global and static user variables */
#define B B_SK
 double B = 0.26;
#define alphar alphar_SK
 double alphar = 7.5;
#define ca0 ca0_SK
 double ca0 = 7e-05;
#define stau stau_SK
 double stau = 10;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ca0_SK", "mM",
 "alphar_SK", "/ms",
 "stau_SK", "ms",
 "B_SK", "mM-cm2/mA-ms",
 "gskbar_SK", "S/cm2",
 "sstau_SK", "ms",
 "cal_SK", "mM",
 "ca_i_SK", "mM",
 "gsk_SK", "S/cm2",
 0,0
};
 static double ca_i0 = 0;
 static double delta_t = 0.01;
 static double q0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ca0_SK", &ca0_SK,
 "alphar_SK", &alphar_SK,
 "stau_SK", &stau_SK,
 "B_SK", &B_SK,
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
 
#define _cvode_ieq _ppvar[5]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"SK",
 "gskbar_SK",
 "sstau_SK",
 "cal_SK",
 "cascale_SK",
 0,
 "gsk_SK",
 0,
 "ca_i_SK",
 "q_SK",
 0,
 0};
 extern Node* nrn_alloc_node_;
 static Symbol* _ca_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gskbar = 0.01;
 	sstau = 20;
 	cal = 0;
 	cascale = 0.1;
 	_prop->param = _p;
 	_prop->param_size = 14;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 6, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 	_ppvar[4]._pval = &nrn_alloc_node_->_area; /* diam */
 prop_ion = need_memb(_ca_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ica */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[1]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _SK_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	ion_reg("k", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 14, 6);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "cvodeieq");
  hoc_register_dparam_semantics(_mechtype, 4, "area");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 SK /Users/eliefarah/ascent/submit/MOD_Files/SK.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   Dca_i = - B * cascale * ica - ( ca_i - ca0 ) / sstau ;
   Dq = alphaq ( _threadargscomma_ ca_i * 1e3 ) * ( 1.0 - q ) - betaq ( _threadargscomma_ ca_i * 1e3 ) * q ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 Dca_i = Dca_i  / (1. - dt*( ( - ( ( 1.0 ) ) / sstau ) )) ;
 Dq = Dq  / (1. - dt*( ( alphaq ( _threadargscomma_ ca_i * 1e3 ) )*( ( ( - 1.0 ) ) ) - ( betaq ( _threadargscomma_ ca_i * 1e3 ) )*( 1.0 ) )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
    ca_i = ca_i + (1. - exp(dt*(( - ( ( 1.0 ) ) / sstau ))))*(- ( ( ( - B )*( cascale ) )*( ica ) - ( ( ( - ca0 ) ) ) / sstau ) / ( ( - ( ( 1.0 ) ) / sstau ) ) - ca_i) ;
    q = q + (1. - exp(dt*(( alphaq ( _threadargscomma_ ca_i * 1e3 ) )*( ( ( - 1.0 ) ) ) - ( betaq ( _threadargscomma_ ca_i * 1e3 ) )*( 1.0 ))))*(- ( ( alphaq ( _threadargscomma_ ca_i * 1e3 ) )*( ( 1.0 ) ) ) / ( ( alphaq ( _threadargscomma_ ca_i * 1e3 ) )*( ( ( - 1.0 ) ) ) - ( betaq ( _threadargscomma_ ca_i * 1e3 ) )*( 1.0 ) ) - q) ;
   }
  return 0;
}
 
double exp1 ( _threadargsprotocomma_ double _lA , double _ld , double _lk , double _lx ) {
   double _lexp1;
  _lexp1 = _lA / exp ( ( 12.0 * log10 ( _lx ) + _ld ) / _lk ) ;
    
return _lexp1;
 }
 
static void _hoc_exp1(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  exp1 ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
double alphaq ( _threadargsprotocomma_ double _lx ) {
   double _lalphaq;
 _lalphaq = exp1 ( _threadargscomma_ 0.00246 , 28.48 , - 4.5 , pow( _lx , 3.0 ) ) ;
   
return _lalphaq;
 }
 
static void _hoc_alphaq(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  alphaq ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double betaq ( _threadargsprotocomma_ double _lx ) {
   double _lbetaq;
 _lbetaq = exp1 ( _threadargscomma_ 0.006 , 60.4 , 35.0 , pow( _lx , 3.0 ) ) ;
   
return _lbetaq;
 }
 
static void _hoc_betaq(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  betaq ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ica = _ion_ica;
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
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
  ica = _ion_ica;
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  ca_i = ca_i0;
  q = q0;
 {
   ca_i = ca0 ;
   q = alphaq ( _threadargscomma_ ca_i * 1e3 ) / ( alphaq ( _threadargscomma_ ca_i * 1e3 ) + betaq ( _threadargscomma_ ca_i * 1e3 ) ) ;
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
  ica = _ion_ica;
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gsk = gskbar * q * q ;
   ik = gsk * ( v - ek ) ;
   cal = ca_i * 1e3 ;
   }
 _current += ik;

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
  ica = _ion_ica;
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
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
  ica = _ion_ica;
  ek = _ion_ek;
 {   state(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(ca_i) - _p;  _dlist1[0] = &(Dca_i) - _p;
 _slist1[1] = &(q) - _p;  _dlist1[1] = &(Dq) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/eliefarah/ascent/submit/MOD_Files/SK.mod";
static const char* nmodl_file_text = 
  "\n"
  ":  SK model based on Aradi and Holmes 1999\n"
  ":\n"
  ":\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX SK\n"
  "	USEION ca READ ica\n"
  "	USEION k READ ek WRITE ik\n"
  "	RANGE gskbar, cal, gsk, sstau, cascale\n"
  "	GLOBAL ca0, stau\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)\n"
  "	(mV) = (millivolt)\n"
  "	(mA) = (milliamp)\n"
  "	(S) = (siemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gskbar = .01	(S/cm2)	: maximum permeability\n"
  "	ca0 = .00007	(mM)\n"
  "	sstau = 20		(ms)\n"
  "	alphar = 7.5	(/ms)\n"
  "	stau = 10		(ms)\n"
  "	cal = 0 (mM)\n"
  "	cascale = .1     ()\n"
  "	B = .26 (mM-cm2/mA-ms)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v		(mV)\n"
  "	ek		(mV)\n"
  "	ik		(mA/cm2)\n"
  "	ica		(mA/cm2)\n"
  "	area		(microm2)\n"
  "      gsk		(S/cm2)\n"
  "}\n"
  "\n"
  "STATE { ca_i (mM) q }\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	gsk = gskbar*q*q\n"
  "	ik = gsk*(v - ek)\n"
  "	cal = ca_i*1e3\n"
  "}\n"
  "\n"
  "DERIVATIVE state {	: exact when v held constant; integrates over dt step\n"
  "	ca_i' = -B*cascale*ica-(ca_i-ca0)/sstau\n"
  "	q' = alphaq(ca_i*1e3)*(1-q)-betaq(ca_i*1e3)*q\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	ca_i = ca0\n"
  "	q = alphaq(ca_i*1e3)/(alphaq(ca_i*1e3)+betaq(ca_i*1e3))\n"
  "}\n"
  "\n"
  "FUNCTION exp1(A (/ms), d, k, x (mM)) (/ms) {\n"
  "	UNITSOFF\n"
  "	exp1 = A/exp((12*log10(x)+d)/k)\n"
  "	UNITSON\n"
  "}\n"
  "\n"
  "FUNCTION alphaq(x (mM)) (/ms) {\n"
  "	alphaq = exp1(0.00246,28.48,-4.5,x^3)\n"
  "}\n"
  "\n"
  "FUNCTION betaq(x (mM)) (/ms) {\n"
  "	betaq = exp1(0.006,60.4,35,x^3)\n"
  "}\n"
  ;
#endif

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
 
#define nrn_init _nrn_init__pump
#define _nrn_initial _nrn_initial__pump
#define nrn_cur _nrn_cur__pump
#define _nrn_current _nrn_current__pump
#define nrn_jacob _nrn_jacob__pump
#define nrn_state _nrn_state__pump
#define _net_receive _net_receive__pump 
 
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
#define pumpbar _p[0]
#define napumpleak _p[1]
#define kpumpleak _p[2]
#define n _p[3]
#define inapump _p[4]
#define ikpump _p[5]
#define nai _p[6]
#define ina _p[7]
#define ik _p[8]
#define v _p[9]
#define _g _p[10]
#define _ion_nai	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ik	*_ppvar[3]._pval
#define _ion_dikdv	*_ppvar[4]._pval
 
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
 "setdata_pump", _hoc_setdata,
 0, 0
};
 /* declare global and static user variables */
#define km2 km2_pump
 double km2 = 67.6;
#define km1 km1_pump
 double km1 = 6.7;
#define kout kout_pump
 double kout = 5;
#define km0 km0_pump
 double km0 = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "km0_pump", "mM",
 "kout_pump", "mM",
 "km1_pump", "mM",
 "km2_pump", "mM",
 "pumpbar_pump", "mA/cm2",
 "napumpleak_pump", "mA/cm2",
 "kpumpleak_pump", "mA/cm2",
 "inapump_pump", "mA/cm2",
 "ikpump_pump", "mA/cm2",
 0,0
};
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "km0_pump", &km0_pump,
 "kout_pump", &kout_pump,
 "km1_pump", &km1_pump,
 "km2_pump", &km2_pump,
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
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"pump",
 "pumpbar_pump",
 "napumpleak_pump",
 "kpumpleak_pump",
 "n_pump",
 0,
 "inapump_pump",
 "ikpump_pump",
 0,
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 11, _prop);
 	/*initialize range parameters*/
 	pumpbar = 0.04;
 	napumpleak = 0;
 	kpumpleak = 0;
 	n = 1.5;
 	_prop->param = _p;
 	_prop->param_size = 11;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* nai */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 	_ppvar[3]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[4]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _fpump_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 11, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "k_ion");
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 pump /Users/eliefarah/ascent/submit/MOD_Files/fpump.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "PUMP";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
 {
   inapump = pumpbar * ( 1.0 / ( 1.0 + pow ( km1 / nai , n ) ) ) ;
   ina = 3.0 * inapump ;
   ik = - 2.0 * inapump ;
   ikpump = ik ;
   kpumpleak = ik ;
   napumpleak = ina ;
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
  nai = _ion_nai;
 initmodel(_p, _ppvar, _thread, _nt);
  }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   inapump = pumpbar * ( 1.0 / ( 1.0 + pow ( km1 / nai , n ) ) ) ;
   ina = ( 3.0 * inapump ) - napumpleak ;
   ik = ( - 2.0 * inapump ) - kpumpleak ;
   ikpump = ik ;
   }
 _current += ina;
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
  nai = _ion_nai;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
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

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/eliefarah/ascent/submit/MOD_Files/fpump.mod";
static const char* nmodl_file_text = 
  "\n"
  "TITLE PUMP\n"
  ": Sodium potassium pump -- Canavier, 1999\n"
  ":\n"
  "\n"
  "UNITS {\n"
  "	(molar) = (1/liter)\n"
  "        (pA) = (picoamp)\n"
  "	(mV) =	(millivolt)\n"
  "        (uS) = (micromho)\n"
  "	(mA) =	(milliamp)\n"
  "	(mM) =	(millimolar)\n"
  "}\n"
  "\n"
  "INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}\n"
  "\n"
  "\n"
  "NEURON {\n"
  "   SUFFIX pump\n"
  "   USEION na READ nai WRITE ina\n"
  "   USEION k WRITE ik\n"
  "   RANGE pumpbar, km, n, inapump, ikpump, kpumpleak, napumpleak : electroneutral sodium accumulation\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "PARAMETER {\n"
  "	  dt (ms)\n"
  "    km0 = 1   (mM)\n"
  "    kout = 5    (mM)\n"
  "    km1 = 6.7   (mM)\n"
  "    km2 = 67.6	(mM)\n"
  "   pumpbar = 0.04   (mA/cm2)\n"
  "   napumpleak = 0 (mA/cm2)\n"
  "   kpumpleak = 0 (mA/cm2)\n"
  "   nai (mM)\n"
  "   n = 1.5\n"
  "   celsius = 35 (degC)\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED{\n"
  "	ina 	(mA/cm2)\n"
  "	ik 	(mA/cm2)\n"
  "	inapump  (mA/cm2)\n"
  "	ikpump (mA/cm2)\n"
  "}\n"
  "\n"
  "\n"
  "INITIAL {\n"
  "\n"
  "inapump = pumpbar*(1/(1+pow(km1/nai,n)))\n"
  "ina = 3.0*inapump\n"
  "ik = -2.0*inapump\n"
  "ikpump = ik\n"
  "\n"
  "kpumpleak = ik\n"
  "napumpleak = ina\n"
  "\n"
  "\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "\n"
  "inapump = pumpbar*(1/(1+pow(km1/nai,n)))\n"
  "ina = (3.0*inapump) - napumpleak\n"
  "ik = (-2.0*inapump) - kpumpleak\n"
  "ikpump = ik\n"
  "\n"
  "}\n"
  ;
#endif

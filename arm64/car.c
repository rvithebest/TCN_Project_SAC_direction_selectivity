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
 
#define nrn_init _nrn_init__car
#define _nrn_initial _nrn_initial__car
#define nrn_cur _nrn_cur__car
#define _nrn_current _nrn_current__car
#define nrn_jacob _nrn_jacob__car
#define nrn_state _nrn_state__car
#define _net_receive _net_receive__car 
#define kin kin__car 
 
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
#define gmax _p[0]
#define gmax_columnindex 0
#define minf _p[1]
#define minf_columnindex 1
#define hinf _p[2]
#define hinf_columnindex 2
#define taum _p[3]
#define taum_columnindex 3
#define tauh _p[4]
#define tauh_columnindex 4
#define mO _p[5]
#define mO_columnindex 5
#define mC _p[6]
#define mC_columnindex 6
#define hO _p[7]
#define hO_columnindex 7
#define hC _p[8]
#define hC_columnindex 8
#define DmO _p[9]
#define DmO_columnindex 9
#define DmC _p[10]
#define DmC_columnindex 10
#define DhO _p[11]
#define DhO_columnindex 11
#define DhC _p[12]
#define DhC_columnindex 12
#define cai _p[13]
#define cai_columnindex 13
#define cao _p[14]
#define cao_columnindex 14
#define ica _p[15]
#define ica_columnindex 15
#define v _p[16]
#define v_columnindex 16
#define _g _p[17]
#define _g_columnindex 17
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
 
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
 static void _hoc_ghkg(void);
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
 "setdata_car", _hoc_setdata,
 "ghkg_car", _hoc_ghkg,
 0, 0
};
#define ghkg ghkg_car
 extern double ghkg( _threadargsprotocomma_ double , double , double , double );
 /* declare global and static user variables */
#define q10 q10_car
 double q10 = 2;
#define taum_exp taum_exp_car
 double taum_exp = 0.92;
#define z z_car
 double z = 2;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gmax_car", 0, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "taum_exp_car", "ms",
 "gmax_car", "S/cm2",
 "taum_car", "ms",
 "tauh_car", "ms",
 0,0
};
 static double delta_t = 0.01;
 static double hC0 = 0;
 static double hO0 = 0;
 static double mC0 = 0;
 static double mO0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "q10_car", &q10_car,
 "taum_exp_car", &taum_exp_car,
 "z_car", &z_car,
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
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"car",
 "gmax_car",
 0,
 "minf_car",
 "hinf_car",
 "taum_car",
 "tauh_car",
 0,
 "mO_car",
 "mC_car",
 "hO_car",
 "hC_car",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 18, _prop);
 	/*initialize range parameters*/
 	gmax = 0.003;
 	_prop->param = _p;
 	_prop->param_size = 18;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
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

 void _car_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
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
  hoc_register_prop_size(_mechtype, 18, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 car /Users/rajuiyer/Desktop/TCN_Project/car.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define FARADAY _nrnunit_FARADAY[_nrnunit_use_legacy_]
static double _nrnunit_FARADAY[2] = {0x1.78e555060882cp+16, 96485.3}; /* 96485.3321233100141 */
 
#define R _nrnunit_R[_nrnunit_use_legacy_]
static double _nrnunit_R[2] = {0x1.0a1013e8990bep+3, 8.3145}; /* 8.3144626181532395 */
static int _reset;
static char *modelname = "Ca R-type channel with medium threshold for activation";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 extern double *_nrn_thread_getelm(SparseObj*, int, int);
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  1
 static int _spth1 = 1;
 static int _cvspth1 = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[4], _dlist1[4]; static double *_temp1;
 static int kin();
 
static int kin (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=2;_i<4;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 minf = 1.0 / ( 1.0 + exp ( - ( v - 3.0 ) / 8.3 ) ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( v + 39.0 ) / 9.2 ) ) ;
   /* ~ mC <-> mO ( minf / taum , ( 1.0 - minf ) / taum )*/
 f_flux =  minf / taum * mC ;
 b_flux =  ( 1.0 - minf ) / taum * mO ;
 _RHS1( 3) -= (f_flux - b_flux);
 
 _term =  minf / taum ;
 _MATELM1( 3 ,3)  += _term;
 _term =  ( 1.0 - minf ) / taum ;
 _MATELM1( 3 ,1)  -= _term;
 /*REACTION*/
  /* ~ hC <-> hO ( hinf / tauh , ( 1.0 - hinf ) / tauh )*/
 f_flux =  hinf / tauh * hC ;
 b_flux =  ( 1.0 - hinf ) / tauh * hO ;
 _RHS1( 2) -= (f_flux - b_flux);
 
 _term =  hinf / tauh ;
 _MATELM1( 2 ,2)  += _term;
 _term =  ( 1.0 - hinf ) / tauh ;
 _MATELM1( 2 ,0)  -= _term;
 /*REACTION*/
   /* mC + mO = 1.0 */
 _RHS1(1) =  1.0;
 _MATELM1(1, 1) = 1;
 _RHS1(1) -= mO ;
 _MATELM1(1, 3) = 1;
 _RHS1(1) -= mC ;
 /*CONSERVATION*/
  /* hC + hO = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= hO ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= hC ;
 /*CONSERVATION*/
   } return _reset;
 }
 
double ghkg ( _threadargsprotocomma_ double _lv , double _lci , double _lco , double _lz ) {
   double _lghkg;
 double _lxi , _lf , _lexi , _lfxi ;
 _lf = R * ( celsius + 273.15 ) / ( _lz * ( 1e-3 ) * FARADAY ) ;
   _lxi = _lv / _lf ;
   _lexi = exp ( _lxi ) ;
   if ( fabs ( _lxi ) < 1e-4 ) {
     _lfxi = 1.0 - _lxi / 2.0 ;
     }
   else {
     _lfxi = _lxi / ( _lexi - 1.0 ) ;
     }
   _lghkg = _lf * ( ( _lci / _lco ) * _lexi - 1.0 ) * _lfxi ;
   
return _lghkg;
 }
 
static void _hoc_ghkg(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  ghkg ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<4;_i++) _p[_dlist1[_i]] = 0.0;}
 minf = 1.0 / ( 1.0 + exp ( - ( v - 3.0 ) / 8.3 ) ) ;
 hinf = 1.0 / ( 1.0 + exp ( ( v + 39.0 ) / 9.2 ) ) ;
 /* ~ mC <-> mO ( minf / taum , ( 1.0 - minf ) / taum )*/
 f_flux =  minf / taum * mC ;
 b_flux =  ( 1.0 - minf ) / taum * mO ;
 DmC -= (f_flux - b_flux);
 DmO += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ hC <-> hO ( hinf / tauh , ( 1.0 - hinf ) / tauh )*/
 f_flux =  hinf / tauh * hC ;
 b_flux =  ( 1.0 - hinf ) / tauh * hO ;
 DhC -= (f_flux - b_flux);
 DhO += (f_flux - b_flux);
 
 /*REACTION*/
   /* mC + mO = 1.0 */
 /*CONSERVATION*/
  /* hC + hO = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<4;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 minf = 1.0 / ( 1.0 + exp ( - ( v - 3.0 ) / 8.3 ) ) ;
 hinf = 1.0 / ( 1.0 + exp ( ( v + 39.0 ) / 9.2 ) ) ;
 /* ~ mC <-> mO ( minf / taum , ( 1.0 - minf ) / taum )*/
 _term =  minf / taum ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 1 ,3)  -= _term;
 _term =  ( 1.0 - minf ) / taum ;
 _MATELM1( 3 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ hC <-> hO ( hinf / tauh , ( 1.0 - hinf ) / tauh )*/
 _term =  hinf / tauh ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 0 ,2)  -= _term;
 _term =  ( 1.0 - hinf ) / tauh ;
 _MATELM1( 2 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
   /* mC + mO = 1.0 */
 /*CONSERVATION*/
  /* hC + hO = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 4;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 4; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 4, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
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
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(_thread[_cvspth1]._pvoid);
   _nrn_destroy_sparseobj_thread(_thread[_spth1]._pvoid);
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  hC = hC0;
  hO = hO0;
  mC = mC0;
  mO = mO0;
 {
   taum = pow( q10 , ( - ( celsius - 22.0 ) / 10.0 ) ) * taum_exp ;
   tauh = pow( q10 , ( - ( celsius - 22.0 ) / 10.0 ) ) * 53.0 ;
    _ss_sparse_thread(&_thread[_spth1]._pvoid, 4, _slist1, _dlist1, _p, &t, dt, kin, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 4; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 ica = gmax * mO * mO * mO * hO * ghkg ( _threadargscomma_ v , cai , cao , z ) ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ica = gmax * mO * mO * mO * hO * ghkg ( _threadargscomma_ v , cai , cao , z ) ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
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
  cai = _ion_cai;
  cao = _ion_cao;
 {  sparse_thread(&_thread[_spth1]._pvoid, 4, _slist1, _dlist1, _p, &t, dt, kin, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 4; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = hO_columnindex;  _dlist1[0] = DhO_columnindex;
 _slist1[1] = mO_columnindex;  _dlist1[1] = DmO_columnindex;
 _slist1[2] = hC_columnindex;  _dlist1[2] = DhC_columnindex;
 _slist1[3] = mC_columnindex;  _dlist1[3] = DmC_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/rajuiyer/Desktop/TCN_Project/car.mod";
static const char* nmodl_file_text = 
  "TITLE Ca R-type channel with medium threshold for activation\n"
  "\n"
  "\n"
  "COMMENT \n"
  "  Kinetics taken from Jeffrey C Magee and Daniel Johnston (1995)\n"
  "  \"Characterization of single voltage-gated Na+ and Ca2+ channels in\n"
  "  apical dendrites of rat CA1 pyramidal neurons\", J. Physiol. 487(1):\n"
  "  67-90.\n"
  "\n"
  "Sterratt et al. 2012\n"
  "http://senselab.med.yale.edu/modeldb/ShowModel.asp?model=144490\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX car\n"
  "    USEION ca READ cai, cao WRITE ica\n"
  "    RANGE gmax, m, h\n"
  "    RANGE minf, hinf, taum, tauh\n"
  "    GLOBAL q10, taum_exp, z\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    (molar) = (/liter)\n"
  "    (mA) = (milliamp)\n"
  "    (nA) = (nanoamp)\n"
  "    (mV) = (millivolt)\n"
  "    (mM) = (millimolar)\n"
  "    (S) = (siemens)\n"
  "    (uS) = (microsiemens)\n"
  "    FARADAY = (faraday) (coulomb)\n"
  "    R = (k-mole) (joule/degC)\n"
  "}\n"
  "\n"
  "PARAMETER {   \n"
  "    gmax = 0.003      (S/cm2) <0,1e9> \n"
  "    : q10  = 3\n"
  "    q10=2  \n"
  "    taum_exp = 0.92  (ms)            : experimentally-measured taum\n"
  "    z = 2                         : valency of Ca ions\n"
  "}  \n"
  "\n"
  "STATE {	mO mC hO hC }    \n"
  "\n"
  "ASSIGNED {               : parameters needed to solve DE\n"
  "    v       (mV)\n"
  "    celsius (degC)\n"
  "    cai     (mM)\n"
  "    cao     (mM)\n"
  "	ica     (mA/cm2)\n"
  "    minf\n"
  "    hinf\n"
  "	taum    (ms)\n"
  "    tauh    (ms)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE kin METHOD sparse\n"
  "	ica = gmax*mO*mO*mO*hO*ghkg(v,cai,cao,z)\n"
  "}\n"
  "\n"
  "INITIAL { \n"
  "    taum = q10^(-(celsius-22(degC))/10(degC))*taum_exp\n"
  "    tauh = q10^(-(celsius-22(degC))/10(degC))*53(ms)\n"
  "    SOLVE kin STEADYSTATE sparse    \n"
  "    ica = gmax*mO*mO*mO*hO*ghkg(v,cai,cao,z)\n"
  "}\n"
  "\n"
  "KINETIC kin {\n"
  "    minf = 1/(1+exp(-(v- 3(mV))/8.3(mV)))\n"
  "    hinf = 1/(1+exp( (v+39(mV))/9.2(mV)))\n"
  "    ~ mC <-> mO (minf/taum, (1-minf)/taum)\n"
  "    ~ hC <-> hO (hinf/tauh, (1-hinf)/tauh)\n"
  "    CONSERVE mC + mO = 1\n"
  "    CONSERVE hC + hO = 1\n"
  "}\n"
  "\n"
  "FUNCTION ghkg(v(mV), ci(mM), co(mM), z) (mV) {\n"
  "    LOCAL xi, f, exi, fxi\n"
  "    f = R*(celsius+273.15)/(z*(1e-3)*FARADAY)\n"
  "    xi = v/f\n"
  "    exi = exp(xi)\n"
  "    if (fabs(xi) < 1e-4) {\n"
  "        fxi = 1 - xi/2\n"
  "    }else{\n"
  "        fxi = xi/(exi - 1)\n"
  "    }\n"
  "    ghkg = f*((ci/co)*exi - 1)*fxi\n"
  "}\n"
  ;
#endif

TITLE minimal model of NMDA receptors

COMMENT
-----------------------------------------------------------------------------

	Minimal kinetic model for glutamate NMDA receptors
	==================================================

  Model of Destexhe, Mainen & Sejnowski, 1994:

	(closed) + T <-> (open)

  The simplest kinetics are considered for the binding of transmitter (T)
  to open postsynaptic receptors.   The corresponding equations are in
  similar form as the Hodgkin-Huxley model:

	dr/dt = alpha * [T] * (1-r) - beta * r

	I = gmax * [open] * B(V) * (V-Erev)

  where [T] is the transmitter concentration and r is the fraction of 
  receptors in the open form.  B(V) represents the voltage-dependent 
  Mg++ block (see Jahr & Stevens J Neurosci 10: 1830, 1990; Jahr & Stevens
  J Neurosci 10: 3178, 1990).

-----------------------------------------------------------------------------

  Based on voltage-clamp recordings of NMDA receptor-mediated currents in rat
  hippocampal slices (Hessler et al., Nature 366: 569-572, 1993), this model 
  was fit directly to experimental recordings in order to obtain the optimal
  values for the parameters (see Destexhe, Mainen and Sejnowski, 1996).

-----------------------------------------------------------------------------

  This mod file includes a mechanism to describe the time course of transmitter
  on the receptors.  The time course is approximated here as a brief pulse
  triggered when the presynaptic compartment produces an action potential.
  The trigger is sensed by the NET_RECEIVE function attached at the
  last, which sets the value of C to Cmax once the trigger is received.

-----------------------------------------------------------------------------

  See details in:

  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for
  computing synaptic conductances based on a kinetic model of receptor binding
  Neural Computation 6: 10-14, 1994.  

  Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of 
  synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition; 
  edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1998, pp. 1-28.

-----------------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDA
	RANGE C, g, gmax, B, lastrelease, TRise, tau
	NONSPECIFIC_CURRENT i
	RANGE Cmax, Cdur, Alpha, Beta, Erev, mg, Deadtime
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	TRise (ms)
	tau   (ms)
	Cmax	= 1	(mM)		: max transmitter concentration
	Erev	= 0	(mV)		: reversal potential
	Prethresh = 0 			: voltage level nec for release
	Deadtime = 1	(ms)		: mimimum time between release events
	gmax		(umho)		: max conductance (100 pS single syn)
	mg	= 1    (mM)		: external magnesium concentration
}

ASSIGNED {
	Alpha	(/ms mM)	: forward (binding) rate
	Beta	(/ms)		: backward (unbinding) rate
	Cdur	(ms)		: transmitter duration (rising phase)
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	C		(mM)		: transmitter concentration
	lastrelease	(ms)		: time of last spike
	B				: magnesium block
}

STATE 
{
	R				: fraction of open channels
}	
	

INITIAL {
	R = 0
	C = 0
	lastrelease = -1000
	Cdur=TRise
	Beta=1/tau
	Alpha=1/Cdur - Beta
}

BREAKPOINT {
	SOLVE states METHOD euler
	B = mgblock(v)		: B is the block by magnesium at this voltage
    g = (gmax * R * B * (Alpha+Beta)) / (Alpha*(1-1/exp(1)))

	i = g*(v - Erev)
}

DERIVATIVE states
{
	evaluateC()
	R'=Alpha * C * (1-R) - Beta * R
}	

	
PROCEDURE evaluateC() 
{ 
	LOCAL q
	q = ((t - lastrelease) - Cdur)		: time since last release ended
	if (q >= 0 && q <= Deadtime && C == Cmax) { : in dead time after release
		C = 0.
	}
}

FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	: from Jahr & Stevens

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

NET_RECEIVE (weight (umho)){

	LOCAL q
	q = ((t - lastrelease) - Cdur)		: time since last release ended

: Spike has arrived, ready for another release?

	if (q > Deadtime) {
		C = Cmax			: start new release
		lastrelease = t
	} 
}	

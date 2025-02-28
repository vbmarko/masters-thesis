
TITLE AMPA Synaptic current


: Declare public variables
NEURON {
	POINT_PROCESS ExpAMPAR
	RANGE  e, i, tau1, tau2, g, mg
	NONSPECIFIC_CURRENT i, iampa
	GLOBAL total
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
	tau1 = 0.5 (ms)			  : tau of alpha synapse for AMPA-R
	tau2 = 1.5 (ms)
	e = 0.0	(mV)
	mg = 1 (mM) 			  : external magnesium concentration
	sh = 0 (mV)

}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	iampa (nA)
	factor
	total (uS)
}

STATE {

	Aampa (uS)
	Gampa (uS)
}

INITIAL {
	LOCAL tp
	total = 0
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	Aampa = 0
	Gampa = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor

}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = Gampa - Aampa
	iampa = g*(v - e)
	i = iampa
}

DERIVATIVE state {
	Aampa' = - Aampa/tau1
	Gampa' = - Gampa/tau2 : - Gampa/tau1	Aampa -> Gampa -> disappear with rate const of 1/tau1.
}

NET_RECEIVE(weight (uS)) {
	state_discontinuity(Aampa, Aampa + weight*factor)
	state_discontinuity(Gampa, Gampa + weight*factor)
	total = total+weight
}


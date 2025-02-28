TITLE A simple passive channel
: This mechanism implements a passive leak current.
: Author: Your Name, Date: YYYY-MM-DD

NEURON {
    SUFFIX leak
    NONSPECIFIC_CURRENT i
    RANGE g_pas, e_pas
}

PARAMETER {
    g_pas = 0.001 (S/cm2)
    e_pas = -65 (mV)
}

ASSIGNED {
    v   (mV)
    i   (mA/cm2)
}

BREAKPOINT {
    i = g_pas*(v - e_pas)
}


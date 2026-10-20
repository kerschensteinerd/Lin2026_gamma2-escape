TITLE Delayed-rectifier Kv current (VG3 amacrine cell) -- sustained component, repositioned to this model's physiological voltage range

COMMENT
Single-gate (n^4), non-inactivating delayed-rectifier Kv current:
ik = gkbar*n^4*(v-ek).

The sustained, 4-AP-sensitive outward component is motivated by Grimes,
Seal, Oesch, Edwards, Diamond (2011), Vis. Neurosci. 28:381-392, Fig. 3E.
The released model sets gkbar from Python (src/vg3_base_functions.py) to
0.02971 S/cm2; calibrate.py checks the leak-subtracted somatic clamp
current after a -90 to 0 mV step against a 993 pA target.

The activation midpoint (-24 mV), slope (3 mV), time constant (2 ms), and
potassium reversal (-85 mV) are modeling assumptions. The activation
curve is positioned within the model's stimulus-driven voltage range;
it is not a direct fit to the voltage dependence in Grimes Fig. 3E.
Grimes did not report a time constant for this current. These assumptions
should be distinguished from the current-magnitude calibration.
ENDCOMMENT

NEURON {
    SUFFIX khh_vg3
    USEION k READ ek WRITE ik
    RANGE gkbar, vhalf, vslope, tau_n, ik
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
}

PARAMETER {
    gkbar  = 0.0 (S/cm2) <0,1e9>
    vhalf  = -24 (mV)
    vslope = 3   (mV)
    tau_n  = 2   (ms)
}

STATE { n }

ASSIGNED {
    v    (mV)
    ek   (mV)
    ik   (mA/cm2)
    ninf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gkbar * n*n*n*n * (v - ek)
}

INITIAL {
    ninf = 1 / (1 + exp(-(v - vhalf) / vslope))
    n = ninf
}

DERIVATIVE states {
    ninf = 1 / (1 + exp(-(v - vhalf) / vslope))
    n' = (ninf - n) / tau_n
}

"""Shared v3 conductance protocol; explicit stimulus and residual inhibition."""
import numpy as np
from neuron import h
from run_fig3_pilot import simulate


def simulate_v3(model,exc_nodes,exc_weights,inh_weights,record_nodes,parameters,
                gaba_remaining,contrast,dt=.1,sample_dt=2.,tstop=3052.):
    """All weights before contrast scaling; fixed same glycine in both groups.

    For the sustained condition, excitation switches reversal from rest to Eexc
    at 50 ms for 3000 ms. Because voltage is initially exactly rest, the dormant
    pre-step shunt produces no baseline current. The active stimulus is a fixed
    excitatory conductance; its amplitude alone scales with contrast. Inhibition
    is a fixed pure shunt at rest throughout, matching the spatial probe.
    """
    tonic=parameters.get('inhibition_protocol','phasic')=='tonic_pure_shunt'
    if parameters.get('excitation_protocol','single_event')!='sustained_3s':
        return simulate(model,exc_nodes,np.asarray(exc_weights)*contrast,
            inh_weights if tonic else np.asarray(inh_weights)*contrast,record_nodes,
            gaba_fraction=1-parameters['glycine_fraction'],gaba_remaining=gaba_remaining,
            tonic=tonic,e_inh=model['rest_v'] if tonic else -60.,dt=dt,sample_dt=sample_dt,tstop=tstop)
    if not tonic: raise ValueError('Sustained v3 stimulus requires explicit fixed pure-shunt inhibition.')
    clamps=[];rest=model['rest_v']
    for index,weight in zip(exc_nodes,np.asarray(exc_weights)*contrast):
        if weight<=0: continue
        clamp=h.SEClamp(model['segments'][int(index)])
        clamp.dur1=50.;clamp.amp1=rest
        clamp.dur2=3000.;clamp.amp2=model['p']['exc_rev']
        clamp.dur3=tstop;clamp.amp3=rest;clamp.rs=1./float(weight)
        clamps.append(clamp)
    for fraction in [(1-parameters['glycine_fraction'])*gaba_remaining,parameters['glycine_fraction']]:
        for index,total in zip(model['inh_nodes'],inh_weights):
            weight=float(total*fraction)
            if weight<=0:continue
            clamp=h.SEClamp(model['segments'][int(index)])
            clamp.dur1=tstop+1.;clamp.amp1=rest;clamp.rs=1./weight
            clamps.append(clamp)
    tv=h.Vector().record(h._ref_t,sample_dt)
    vv=[h.Vector().record(model['segments'][int(i)]._ref_v,sample_dt) for i in record_nodes]
    h.CVode().active(0);h.dt=dt;h.finitialize(rest);h.continuerun(tstop)
    return np.asarray(tv).copy(),np.array([np.asarray(v) for v in vv])

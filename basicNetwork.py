# -*- coding: utf-8 -*-

import os
import nest
import pylab
import numpy as np
import datetime
import nest.raster_plot

nest.ResetKernel()

# print dir(nest) # All the possible nest command available
# nest.GetStatus(neuron, "I_e")
# nest.GetStatus(neuron, ["V_reset", "V_th"])
# nest.SetStatus(neuron, {"I_e": 376.0})
# nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})
# nest.SetStatus(neuron2 , {"I_e": 370.0})

tSimu = 1000.0
timeID = str(datetime.datetime.now()).replace('-','_').replace(' ','_').replace(':','_')[:-7]
print timeID
noise_exc = nest.Create("poisson_generator", 10, {"rate": 70.0, "start": 0.00, "stop": tSimu})
noise_inh = nest.Create("poisson_generator", 10, {"rate": 45.0, "start": 0.00, "stop": 500.0})

parrot_exc = nest.Create("parrot_neuron", 10)
parrot_inh = nest.Create("parrot_neuron", 10)

pop_exc = nest.Create("iaf_psc_alpha", 10, params={"I_e": 0.0, "V_m": -70.0, "E_L": -70.0, "C_m": 150.0, "tau_m": 10.0, "t_ref": 0.002,
                                                  "V_th": -55.0, "V_reset": -70.0, "tau_syn_ex": 1.0, "tau_syn_in": 1.0, "V_min": -90.0,}) 
pop_inh = nest.Create("iaf_psc_alpha", 10, params={"I_e": 0.0, "V_m": -70.0, "E_L": -70.0, "C_m": 150.0, "tau_m": 10.0, "t_ref": 0.002,
                                                  "V_th": -55.0, "V_reset": -70.0, "tau_syn_ex": 1.0, "tau_syn_in": 1.0, "V_min": -90.0,}) 

pop_out = nest.Create("iaf_psc_alpha", 1, params={"I_e": 200.0,   #  I_e	 double	- Constant external input current in pA.  
                                                   "V_m": -70.0,     #  V_m	 double	- Membrane potential in mV  
                                                   "E_L": -70.0,     #  E_L	 double	- Resting membrane potential in mV.  
                                                   "C_m": 150.0,     #  C_m	 double	- Capacity of the membrane in pF  
                                                   "tau_m": 10.0,    #  tau_m	 double	- Membrane time constant in ms.  
                                                   "t_ref": 0.002,   #  t_ref	 double	- Duration of refractory period in ms.  
                                                   "V_th": -55.0,    #  V_th	 double	- Spike threshold in mV.  
                                                   "V_reset": -70.0, #  V_reset   double	- Reset potential of the membrane in mV.  
                                                   "tau_syn_ex": 1.0,#  tau_syn_ex double	- Rise time of the excitatory synaptic alpha function in ms.  
                                                   "tau_syn_in": 1.0,#  tau_syn_in double	- Rise time of the inhibitory synaptic alpha function in ms.  
                                                   "V_min": -90.0,}) #  V_min	 double	- Absolute lower value for the membrane potential.  

# connection of the inputs
nest.Connect(pre=noise_exc,post=parrot_exc)
nest.Connect(pre=noise_inh,post=parrot_inh)                                
nest.Connect(pre=parrot_exc, post=pop_exc, syn_spec={"weight": 100., "delay":1.0}, conn_spec={'rule':'one_to_one'}, model=None)
nest.Connect(pre=parrot_inh, post=pop_inh, syn_spec={"weight": 100., "delay":1.0},  conn_spec={'rule':'one_to_one'}, model=None)

nest.Connect(pre=pop_exc, post=pop_out, syn_spec={"weight": 75., "delay":1.0}, conn_spec={'rule':'fixed_indegree','indegree': 3}, model=None)
nest.Connect(pre=pop_inh, post=pop_out, syn_spec={"weight": -100., "delay":1.0},  conn_spec={'rule':'fixed_indegree','indegree': 10}, model=None)

multi_out = nest.Create("multimeter", params={"withtime":True, "record_from":["V_m"]})
spike_out = nest.Create("spike_detector", params={"withgid": True, "withtime": True,"label": 'pop_out_'+timeID, "to_file": False, 'start': 0.00 ,'stop': tSimu})
nest.Connect(multi_out, pop_out)
nest.Connect(pop_out, spike_out)

multi_exc = nest.Create("multimeter", params={"withtime":True, "record_from":["V_m"]})
spike_exc = nest.Create("spike_detector", params={"withgid": True, "withtime": True,"label": 'pop_exc_'+timeID, "to_file": False, 'start': 0.00 ,'stop': tSimu})
nest.Connect(multi_exc, pop_exc)
nest.Connect(pop_exc, spike_exc)

multi_inh = nest.Create("multimeter", params={"withtime":True, "record_from":["V_m"]})
spike_inh = nest.Create("spike_detector", params={"withgid": True, "withtime": True,"label": 'pop_inh_'+timeID, "to_file": False, 'start': 0.00 ,'stop': tSimu})
nest.Connect(multi_inh, pop_inh)
nest.Connect(pop_inh, spike_inh)

nest.Simulate(tSimu)



frr_out = nest.GetStatus(spike_out, 'n_events')[0]*1000 / float(tSimu)
frr_exc = nest.GetStatus(spike_exc, 'n_events')[0]*1000 / float(tSimu*10)
frr_inh = nest.GetStatus(spike_inh, 'n_events')[0]*1000 / float(tSimu*10)

print "---------------"
print "FRR out = ", frr_out, "Hz"
print "FRR exc = ", frr_exc, "Hz"
print "FRR inh = ", frr_inh, "Hz"
print "---------------"

dmm_out = nest.GetStatus(multi_out)[0]
Vms_out = dmm_out["events"]["V_m"]
tsm_out = dmm_out["events"]["times"]

dS_out = nest.GetStatus(spike_out,keys="events")[0]
evs_out = dS_out["senders"]
tsd_out = dS_out["times"]

dmm_exc = nest.GetStatus(multi_exc)[0]
Vms_exc = dmm_exc["events"]["V_m"]
tsm_exc = dmm_exc["events"]["times"]
dS_exc = nest.GetStatus(spike_exc,keys="events")[0]
evs_exc = dS_exc["senders"]
tsd_exc = dS_exc["times"]

dmm_inh = nest.GetStatus(multi_inh)[0]
Vms_inh = dmm_inh["events"]["V_m"]
tsm_inh = dmm_inh["events"]["times"]
dS_inh = nest.GetStatus(spike_inh,keys="events")[0]
evs_inh = dS_inh["senders"]
tsd_inh = dS_inh["times"]

pylab.figure(1,figsize=[18, 6])

pylab.subplot(2, 3, 1)
pylab.title('out population')
pylab.plot(tsm_out, Vms_out, '-r')

pylab.subplot(2, 3, 2)
pylab.title('excitatory population')
pylab.plot(tsm_exc, Vms_exc, '-b')

pylab.subplot(2, 3, 3)
pylab.title('inhibitory population')
pylab.plot(tsm_inh, Vms_inh, '-g')

pylab.subplot(2, 3, 4)
pylab.plot(tsd_out, evs_out, "|r")

pylab.subplot(2, 3, 5)
pylab.plot(tsd_exc, evs_exc, "|b")

pylab.subplot(2, 3, 6)
pylab.plot(tsd_inh, evs_inh, "|g")

# pylab.savefig('eZKUFYV')
pylab.show()



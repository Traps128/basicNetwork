# -*- coding: utf-8 -*-

import os
import nest
import pylab
import numpy as np
import nest.raster_plot
# print dir(nest) # All the possible nest command available
# nest.GetStatus(neuron, "I_e")
# nest.GetStatus(neuron, ["V_reset", "V_th"])
# nest.SetStatus(neuron, {"I_e": 376.0})
# nest.SetStatus(multimeter, {"withtime":True, "record_from":["V_m"]})
# nest.SetStatus(neuron2 , {"I_e": 370.0})

#tSimu = 1000.0
#
#pop_out = nest.Create("iaf_psc_alpha", params={"I_e": 370.0,     #  I_e	 double	- Constant external input current in pA.  
#                                                   "V_m": -70.0,     #  V_m	 double	- Membrane potential in mV  
#                                                   "E_L": -70.0,     #  E_L	 double	- Resting membrane potential in mV.  
#                                                   "C_m": 200.0,     #  C_m	 double	- Capacity of the membrane in pF  
#                                                   "tau_m": 10.0,    #  tau_m	 double	- Membrane time constant in ms.  
#                                                   "t_ref": 0.002,   #  t_ref	 double	- Duration of refractory period in ms.  
#                                                   "V_th": -55.0,    #  V_th	 double	- Spike threshold in mV.  
#                                                   "V_reset": -70.0, #  V_reset   double	- Reset potential of the membrane in mV.  
#                                                   "tau_syn_ex": 1.0,#  tau_syn_ex double	- Rise time of the excitatory synaptic alpha function in ms.  
#                                                   "tau_syn_in": 1.0,#  tau_syn_in double	- Rise time of the inhibitory synaptic alpha function in ms.  
#                                                   "V_min": -90.0,}) #  V_min	 double	- Absolute lower value for the membrane potential.  
#                                    
#
#multimeter = nest.Create("multimeter", params={"withtime":True})
#spikedetector = nest.Create("spike_detector", params={"withgid": True, "withtime": True})
#
#nest.Connect(multimeter, pop_out)
#nest.Connect(pop_out, spikedetector)
#
#nest.Simulate(tSimu)
#    
#frr = nest.GetStatus(spikedetector, 'n_events')[0]*1000 / float(tSimu)
#
#print "---------------"
#print "FRR = ", frr, "Hz"
#print "---------------"
#
#print nest.GetStatus(pop_out)

nest.help(nest.raster_plot)
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 19:19:06 2020

@author: rh17872 - Rachel Humphries

Stimulates increasing numbers of AMPA/NMDA inputs on each dendritic section in the SR/SLM region.

Produces data used in Figures 3 and 4 of paper: 

Rachel Humphries, Jack R. Mellor, Cian O'Donnell, 
Acetylcholine Boosts Dendritic NMDA Spikes in a CA3 Pyramidal Neuron Model,
Neuroscience, 2021,ISSN 0306-4522
https://doi.org/10.1016/j.neuroscience.2021.11.014

"""

import numpy as np
from neuron import h
from neuron.units import ms, mV
h.load_file('stdrun.hoc')
import ca3_synapse_functions as csf
import pickle

pickle_file = "nonlin_outputs.pkl" #pickle file to save outputs
regions = ["SR","SLM"] 

'''Synapse parameters'''
syn_ISIs = 1 #interspike interval
syn_density = 1 #ums between each synapse 0,1,2,5
syn_NAratios = [0.5] #AMPA to NMDA ratio
experiment = [str(x)+" AMPA:NMDA"  for x in syn_NAratios] #for testing different AMPA:NMDA ratios

n_ac_syn = [5,16] #NMDA rise and decay times in SR dendrites
a_ac_syn = [0.5,1.5] #AMPA rise and decay times in SR dendrites
n_pp_syn = [5,16] #NMDA rise and decay times in SLM dendrites
a_pp_syn = [0.5,1.5] #AMPA rise and decay times in SLM dendrites
total_cond_ac = 0.0035*0.5
total_cond_pp = 0.0035*0.5

a_ac_weights, n_ac_weights = csf.calculate_synapse_weights(syn_NAratios,total_cond_ac, "AC")
a_pp_weights, n_pp_weights = csf.calculate_synapse_weights(syn_NAratios,total_cond_pp, "PP")

for i,x in enumerate(a_ac_weights) :
    print ("AC AMPA:NMDA ratios :", x/n_ac_weights[i] )
    print ("PP AMPA:NMDA ratios :", a_pp_weights[i]/n_pp_weights[i])
    
'''Stimulation parameters'''
total_syns = 20 # must be an even number
step = 1
num_synapses = range(step,total_syns+1,step)

syn_density = 1 #number of ums between synapses
stim_num = 1 #number of stims
stim_start = 810 #start of stim (ms)
stim_noise = 0 #noise of stim
stim_interval = 20 #interval between stims (20ms = 50hz)

'''Parameters from optimisation'''
optimiser = [5.557118781038023e-11, 0.0011114634012579765, 0.0011580392555217256, 3.472924393529932e-06, 4.403602806593143e-10, 4.616466181142678e-09, 5.583387888166266e-06, 5.5750027000813285, -29.999999524116564, 30.0, 30.0, -48.167858484536445, 0.5, 0.5, 0.5, 0.20000188521088813]

gka =optimiser[0]#random.uniform(0.01,1) #0.8 #0.02
gkm = optimiser[1]#random.uniform(0.001,0.1)#0.1 #0.017 #0.03
gkca =optimiser[2]#random.uniform(0.001,0.1)#0.01 #0.001
gkir =optimiser[3]#random.uniform(0.0001,0.01)#0.001 #0.00015 #1.44e-4
gpas =optimiser[4]#random.uniform(0.0004,0.04)#0.0004 #4e-8
gkdr = optimiser[5]
gih = optimiser[6]

ka_sh = optimiser[7]
km_sh = optimiser[8]#random.uniform(0,25)
ih_sh =  optimiser[9]#random.uniform(0,25)
kir_sh = optimiser[10]#random.uniform(0,25)
epas = optimiser[11] #-60.

ka_ach_block = optimiser[12] #random.uniform(0.,0.2)
km_ach_block = optimiser[13] #random.uniform(0.,0.2)
kca_ach_block = optimiser[14] #random.uniform(0.,0.2)
kir_ach_block = optimiser[15]#random.uniform(0.,0.2)

gkd = 0.0

'''Set up neuron model parameters and conductances'''

h('{Vrest = -75.}')
h('{vrest_val = -75.}')
h('{tstop=1000}')
h('{epas = -70}')
h('{epas_val = -70}')
#h('{Rm = 50740}')
#h('{rm_val = 50740}')
h('{Cm    = 0.7}')
h('{RaAll= 150}')
h('{gpas = 1/25370}')
h('{gpas_val = 1/25370}')
h('AXONM = 5')
h('{gna =  0.0}')
h('{gkdr = 0.005}') #0.00518 #0.125
h('{gkdr_val = 0.005}') #0.00518 #0.125
h('{KMULT =  0.02}') #Ka conductance
h('{gka_val =  0.02}') #Ka conductance
h('{gkm=0.017}')
h('{gkm_val=0.017}')
h('{gkir=1.44e-05}')
h('{gkir_val=1.44e-05}')
h('{gkd=0.0}')
h('{gkd_val=0.0}')
h('{gc=1.e-5}')
h('{gcal=gc}') #1.0659e-5 #0.000507
h('{gcal_val=gc}') #1.0659e-5 #0.000507
h('{gcat=gc}') #5.984e-7 #4.554e-7
h('{gcat_val=gc}') #5.984e-7 #4.554e-7
h('{gcan=gc}') #3.791e-5 #0.000165
h('{gcan_val=gc}') #3.791e-5 #0.000165
h('{gKc=5e-5}') #0.000111 #0.00412
h('{gKc_val=5e-5}')
#h('{gkc_val=5e-5}')
h('{gkcas=0.001}') 
h('{gkcas_val=0.001}')
h('{gahp=0.0001}') #0.000511 #0.00179
h('{gahp_val=0.0001}')
h('{ghd=0.00001}') #6.529e-5 #9.578e-5
h('{ghd_val=0.00001}')
h('km_sh = 0')
h('km_sh_val = 0')
h('ka_sh = 0')
h('ka_sh_val = 0')
h('kir_sh = 0')
h('kir_sh_val = 0')   
h('ih_sh = 0')
h('ih_sh_val = 0')
h('na_block = 1.')
h('na_block_val = 1.')
h('ka_block = 1.')
h('ka_block_val = 1.')

h.load_file('ca3b-cell1zr-fig9b.hoc')

h.gpas_val = gpas
h.epas_val = epas
#h.rm_val = rm

h.gka_val = gka
h.ka_block_val = 1.
if gka == 0.:
    h.ka_block_val = 0.
h.na_block_val = 1.
if h.gna == 0.:
    h.na_block_val = 0.

h.gkm_val = gkm
h.gkcas_val = gkca
h.gkir_val = gkir
h.gkdr_val = gkdr
h.ghd_val = gih

h.km_sh_val = km_sh
h.ka_sh_val= ka_sh
h.kir_sh_val= kir_sh
h.ih_sh_val = ih_sh

h.dt = 1.
dt = 1.

'''Divide dendrites into regions based on distance from soma'''
h.distance(sec=h.soma[0]) #set the origin at the soma
sl_sr_bound = 150
sr_slm_bound = 400
rel = 0.5
sl_dends = []
sr_dends = []
slm_dends = []
sl_dends, sr_dends, slm_dends = csf.divide_dends(rel,sl_sr_bound,sr_slm_bound,sl_dends,sr_dends,slm_dends)
so_dends = range(int(h.numbasal)) #basal dendrites
all_dends = [so_dends, sl_dends, sr_dends, slm_dends]

'''Identify dendrite sections longer than 20um to use in simulation'''
so_dend_200 = []
for d in so_dends :
    if h.dendrite[d].L > 20 :
        so_dend_200.append([d]) 
sr_dend_200 = [[d] for d in sr_dends if h.apical_dendrite[d].L > 20] 
slm_dend_200 = [[d] for d in slm_dends if h.apical_dendrite[d].L > 20]
print ("SR dends: ", sr_dend_200)
print ("SLM dends: ", slm_dend_200)
all_dends_200 = []
all_dends_region = []
if "SO" in regions :
    all_dends_200.append(so_dend_200)
    all_dends_region.append(so_dends)
if "SR" in regions :
    all_dends_200.append(sr_dend_200)
    all_dends_region.append(sr_dends)
if "SLM" in regions :
    all_dends_200.append(slm_dend_200)
    all_dends_region.append(slm_dends)

'''Calculate dendrite diameters'''
dend_diams = csf.get_dend_diameters(all_dends_200, regions)
                
'''Set up lists'''
soma_volts = [[[[] for i in experiment] for j in x]for x in all_dends_200] 
dend_volts = [[[[] for i in experiment] for j in x]for x in all_dends_200]
soma_peaks = [[[[] for i in experiment] for j in x]for x in all_dends_200]
dend_peaks = [[[[] for i in experiment] for j in x]for x in all_dends_200]
dend_distances = [[[] for j in x]for x in all_dends_200]
irs = [[[]for j in x] for x in all_dends_200]
dend_centers = [[[]for j in x]for x in all_dends_200]
nmda_currs = [[[[] for i in experiment] for j in x]for x in all_dends_200]
ampa_soma_volts = [[[[] for i in experiment] for j in x]for x in all_dends_200]
ampa_dend_volts = [[[[] for i in experiment] for j in x]for x in all_dends_200] 
ampa_soma_peaks = [[[[] for i in experiment] for j in x]for x in all_dends_200] 
ampa_dend_peaks = [[[[] for i in experiment] for j in x]for x in all_dends_200] 


'''Run simulation'''        
for r,region in enumerate(all_dends_200) : #loops through region
    print ("Region: ", regions[r])
    for j,dends in enumerate(region) : #loops through dendrite group
        print ("Dendrites: ", dends)
        if regions[r] == "SO" :
            pos_per_dend, mid_dend,mid_syn_pos,ordered_pos = csf.synapse_positioning(h.dendrite,dends,syn_density,num_synapses[-1]) #finds dendrite positions on dendrite section    
        else :
            pos_per_dend, mid_dend,mid_syn_pos,ordered_pos = csf.synapse_positioning(h.apical_dendrite,dends,syn_density,num_synapses[-1]) #finds dendrite positions on dendrite section    
        for s,nar in enumerate(syn_NAratios) : #loops through synapse densities.
            print ("Synapse NA ratio: ", experiment[s])
            for n in num_synapses : #loops through number of synapses to add
                print ("Combined synapses: AMPA+NMDA" , n)
                glu_syns = []
                glu_netstim = []
                glu_netcon = []
                for ns in range(n) : #loops through each synapse to add - adds both AMPA + NMDA synapse
                    d = dends[ordered_pos[ns][0]]
                    p = pos_per_dend[ordered_pos[ns][0]][ordered_pos[ns][1]]
                    #print "Synapses: ", d, p
                    if regions[r] == "SO" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_nmda_baker(h.dendrite[d],p,n_ac_weights[s],n_ac_syn[0],n_ac_syn[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise) #stim_start+(syn_ISIs[0]*ns)
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.dendrite[d],p,a_ac_weights[s],a_ac_syn[0],a_ac_syn[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                    if regions[r] == "SR" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_nmda_baker(h.apical_dendrite[d],p,n_ac_weights[s],n_ac_syn[0],n_ac_syn[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.apical_dendrite[d],p,a_ac_weights[s],a_ac_syn[0],a_ac_syn[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                    if regions[r] == "SLM" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_nmda_baker(h.apical_dendrite[d],p,n_pp_weights[s],n_pp_syn[0],n_pp_syn[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.apical_dendrite[d],p,a_pp_weights[s],a_pp_syn[0],a_pp_syn[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)

                soma_v = h.Vector().record(h.soma[0](0.5)._ref_v,dt) #set up hoc vector to record somatic voltage
                # set up hoc vector to record dendritic voltage
                if regions[r] == "SO" :
                    dend_v = h.Vector().record(h.dendrite[dends[ordered_pos[0][0]]](pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])._ref_v,dt)
                else:
                    dend_v = h.Vector().record(h.apical_dendrite[dends[ordered_pos[0][0]]](pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])._ref_v,dt)
                    
                t = h.Vector().record(h._ref_t,dt) #record time
                nmda_sing_syns = [] # for the individual synaptic currents
                for syn in glu_syns[0::2]: #for recording NMDA currents
                    tmp = h.Vector()
                    tmp.record(syn._ref_i,dt)
                    nmda_sing_syns.append(tmp)
                
                h.fig9b() #calls the function in the hoc file to run the simulation
                soma_v = np.array(soma_v)
                dend_v = np.array(dend_v)
                
                #saves data to lists            
                soma_volts[r][j][s].append(soma_v)
                dend_volts[r][j][s].append(dend_v)
                soma_peaks[r][j][s].append(max(soma_v[int(stim_start/dt):]))
                dend_peaks[r][j][s].append(max(dend_v[int(stim_start/dt):]))
                
                inmda_total = np.zeros(len(t))
                for curr in nmda_sing_syns:
                    inmda_total += curr.as_numpy()
                nmda_currs[r][j][s].append(inmda_total)
            
            '''Repeats simulation with only AMPA synapses'''
            for n in num_synapses : #loops through number of synapses to add
                print ("Only AMPA" , n)
                glu_syns = []
                glu_netstim = []
                glu_netcon = []
                for ns in range(n) : #loops through each synapse to add
                    d = dends[ordered_pos[ns][0]]
                    p = pos_per_dend[ordered_pos[ns][0]][ordered_pos[ns][1]]
                    #print "Synapses: ", d, p
                    if regions[r] == "SO" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.dendrite[d],p,a_ac_weights[s],a_ac_syn[0],a_ac_syn[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                    if regions[r] == "SR" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.apical_dendrite[d],p,a_ac_weights[s],a_ac_syn[0],a_ac_syn[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                    if regions[r] == "SLM" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.apical_dendrite[d],p,a_pp_weights[s],a_pp_syn[0],a_pp_syn[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)

                soma_v = h.Vector().record(h.soma[0](0.5)._ref_v,dt)
                if regions[r] == "SO" :
                    dend_v = h.Vector().record(h.dendrite[dends[ordered_pos[0][0]]](pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])._ref_v,dt)
                else:
                    dend_v = h.Vector().record(h.apical_dendrite[dends[ordered_pos[0][0]]](pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])._ref_v,dt)
                    
                t = h.Vector().record(h._ref_t,dt)
                
                h.fig9b() #calls the function in the hoc file to run the simulation
                soma_v = np.array(soma_v)
                dend_v = np.array(dend_v)
                
                #saves data to lists
                ampa_soma_volts[r][j][s].append(soma_v)
                ampa_dend_volts[r][j][s].append(dend_v)
                ampa_soma_peaks[r][j][s].append(max(soma_v[int(stim_start/dt):]))
                ampa_dend_peaks[r][j][s].append(max(dend_v[int(stim_start/dt):]))
                
        freq=0
        rin= h.Impedance() # does not measure resting Rin here
        rin.compute(freq,1)
        if regions[r] == "SO" :
            input_imp = rin.input(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]],sec=h.dendrite[dends[ordered_pos[0][0]]])        
        else:
            input_imp = rin.input(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]],sec=h.apical_dendrite[dends[ordered_pos[0][0]]])
        irs[r][j].append(input_imp)
        
        if regions[r] == "SO" :
            dend_distances[r][j].append(-h.distance(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]],sec=h.dendrite[dends[ordered_pos[0][0]]]))
        else :
            dend_distances[r][j].append(h.distance(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]],sec=h.apical_dendrite[dends[ordered_pos[0][0]]]))
        
        dend_centers[r][j].append(dends[ordered_pos[0][0]])
        dend_centers[r][j].append(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])
        
for i,r in enumerate(irs):
    for j,d in enumerate(r) :
        print (regions[i], "IRs ", all_dends_200[i][j], d )

time = np.array(t)

with open(pickle_file, 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([all_dends_200,soma_volts, dend_volts, dt, time, irs, soma_peaks, dend_peaks, nmda_currs, ampa_soma_peaks,ampa_dend_peaks,ampa_soma_volts, ampa_dend_volts,dend_distances,sr_dend_200,slm_dend_200,dend_diams,so_dends,sl_dends,sr_dends,slm_dends], f)
f.close()
''' 
Pickle file data

all_dends_200 : dendrite sections used in simulation
soma_volts : somatic voltages of each simulation (mv)
dend_volts : dendritic voltages of stimulated dendrite (mv)
dt : simulation timestep (ms)
time : time
irs : input resistance during simulation (Mohm)
soma_peaks : peak of each somatic voltage trace (mv)
dend_peaks : peak of each dendritic voltage trace (mv)
nmda_currs : NMDA current 
ampa_soma_peaks : peak of AMPA only simulations recorded at soma (mv) 
ampa_dend_peaks : peak of AMPA only simulations recorded in dendrites (mv)
ampa_soma_volts : somatic voltages of each AMPA-only simulation (mv)
ampa_dend_volts : AMPA-only dendritic voltages of stimulated dendrite (mv)
dend_distances : distance of stimulated dendrites from soma (um)
sr_dend_200 : SR dendrite sections used in simulation 
slm_dend_200 : SLM dendrite sections used in simulation   
dend_diams : dendrite diameters 
so_dends : all SO dendrite sections
sl_dends : all SL dendrite sections
sr_dends : all SR dendrite sections
slm_dends : all SLM dendrite sections  
  
'''

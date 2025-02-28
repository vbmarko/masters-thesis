# -*- coding: utf-8 -*-
"""
Created on Fri May 14 15:30:54 2021

@author: rh17872 - Rachel Humphries

Measures resting input resistance of each dendritic section with different potassium channel inhibitions.  

Outputs resting input resistance data used in Figures 4, 5 & 6 of paper: 

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

pickle_file = "Rin_rest.pkl"
experiment = ["Control", "ACh", "Ka" r'$\downarrow$', "Km" r'$\downarrow$',"Kca" r'$\downarrow$',"Kir" r'$\downarrow$']
save_exps = ["Control", "ACh", "Ka", "Km", "Kca", "Kir"]
regions = ["SR", "SLM"]

'''Optimiser parameters'''
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

'''Syanpse parameters''' 
step = 1
total_syns = 20
num_synapses = range(step,total_syns+1,step)
syn_density = 1
syn_ISIs = [1] 
syn_density = 1 #ums between each synapse 0,1,2,5
syn_NAratios = [0.5] #AMPA to NMDA ratio

'''Set up neuron model'''

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

h.gkm_val = gkm
h.gkcas_val = gkca
h.gkir_val = gkir
h.gkdr_val = gkdr
h.ghd_val = gih

h.km_sh_val = km_sh
h.ka_sh_val= ka_sh
h.kir_sh_val= kir_sh
h.ih_sh_val = ih_sh

h.dt = 1000.
dt = 1000.

def setup_potassium_channels(save_exps,pot_channel,pot_conductance,pot_block) :
    if pot_channel in save_exps :
        pot_cond_list= [pot_conductance]*(len(save_exps)) #0.02 / 0.0068 /0.8 (*40)
        #pot_cond_list.append(pot_conductance*pot_block)
        pot_cond_list = [0. if i==save_exps.index(pot_channel) else x for i,x in enumerate(pot_cond_list)]
        pot_cond_list = [pot_conductance*pot_block if save_exps[i]=="ACh" else x for i,x in enumerate(pot_cond_list)]
    else :
        pot_cond_list = [pot_conductance]*(len(save_exps))
        pot_cond_list = [pot_conductance*pot_block if save_exps[i]=="ACh" else x for i,x in enumerate(pot_cond_list)]
        #pot_cond_list.append(pot_conductance*pot_block)
    return pot_cond_list
    
gkas = setup_potassium_channels(save_exps,"Ka",gka,1.) #ka_ach_block has to be set during simulation in the hoc file
gkms = setup_potassium_channels(save_exps,"Km",gkm,km_ach_block)
gkcas = setup_potassium_channels(save_exps,"Kca",gkca,kca_ach_block)
gkirs = setup_potassium_channels(save_exps,"Kir",gkir,kir_ach_block)

print ("Potassium conductances: ", gkas, gkms, gkcas, gkirs)

'''Divide apical dendrites'''
h.distance(sec=h.soma[0]) #set the origin at the soma
sl_dends = []
sr_dends = []
slm_dends = []
sl_sr_bound = 150
sr_slm_bound = 400
rel = 0.5
sl_dends, sr_dends, slm_dends = csf.divide_dends(rel,sl_sr_bound,sr_slm_bound,sl_dends,sr_dends,slm_dends)
so_dends = range(int(h.numbasal)) #basal dendrites
all_dends = [so_dends, sl_dends, sr_dends, slm_dends]

'''Dendrite sections longer than 20um'''
so_dend_200 = []
for d in so_dends : 
    if h.dendrite[d].L > 20:
        so_dend_200.append([d]) 
sr_dend_200 = [[d] for d in sr_dends if h.apical_dendrite[d].L > 20] 
slm_dend_200 = [[d] for d in slm_dends if h.apical_dendrite[d].L > 20]
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

'''Set up lists'''
soma_volts = [[[[] for i in experiment] for j in x]for x in all_dends_200] 
dend_volts = [[[[] for i in experiment] for j in x]for x in all_dends_200]
dend_distances = [[[] for j in x]for x in all_dends_200]
irs = [[[]for j in x] for x in all_dends_200]
soma_irs = [[[]for j in x] for x in all_dends_200]
dend_centers = [[[]for j in x]for x in all_dends_200]

'''Run simulation - no synaptic input'''        
for r,region in enumerate(all_dends_200) : #loops through region
    print ("Region: ", regions[r])
    for j,dends in enumerate(region) : #loops through dendrite group
        print ("Dendrites: ", dends)
        if regions[r] == "SO" :
            pos_per_dend, mid_dend,mid_syn_pos,ordered_pos = csf.synapse_positioning(h.dendrite,dends,syn_density,num_synapses[-1])    
        else :
            pos_per_dend, mid_dend,mid_syn_pos,ordered_pos = csf.synapse_positioning(h.apical_dendrite,dends,syn_density,num_synapses[-1])    
        for s,nar in enumerate(experiment) : #loops through potassium channel blocks.
            print ("Experiment: ", save_exps[s])
            h.gka_val = gkas[s]
            h.ka_block_val = 1.
            if gkas[s] == 0.:
                h.ka_block_val = 0.
            if nar == "ACh":
                h.ka_block_val = ka_ach_block
            h.na_block_val = 1.
            if h.gna == 0.:
                h.na_block_val = 0.
            h.gkm_val = gkms[s]
            h.gkcas_val = gkcas[s]
            h.gkir_val = gkirs[s]
            
            soma_v = h.Vector().record(h.soma[0](0.5)._ref_v,dt)
            if regions[r] == "SO" :
                dend_v = h.Vector().record(h.dendrite[dends[ordered_pos[0][0]]](pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])._ref_v,dt)
            else:
                dend_v = h.Vector().record(h.apical_dendrite[dends[ordered_pos[0][0]]](pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])._ref_v,dt)
                        
            t = h.Vector().record(h._ref_t,dt)
            
            h.fig9b() #calls the function in the hoc file to run the simulation
            
            soma_v = np.array(soma_v)
            dend_v = np.array(dend_v)
        
            soma_volts[r][j][s].append(soma_v)
            dend_volts[r][j][s].append(dend_v)
            
            '''Compute input resistance in each dendrite'''
            freq=0
            rin= h.Impedance()
            rin.compute(freq,1)
            if regions[r] == "SO" :
                input_imp = rin.input(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]],sec=h.dendrite[dends[ordered_pos[0][0]]])        
            else:
                input_imp = rin.input(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]],sec=h.apical_dendrite[dends[ordered_pos[0][0]]])
                soma_ir = rin.input(0.5,sec=h.soma[0])
            
            irs[r][j].append(input_imp) #dendritic IR
            soma_irs[r][j].append(soma_ir)
        
        if regions[r] == "SO" :
            dend_distances[r][j].append(-h.distance(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]],sec=h.dendrite[dends[ordered_pos[0][0]]]))
        else :
            dend_distances[r][j].append(h.distance(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]],sec=h.apical_dendrite[dends[ordered_pos[0][0]]]))
        
        dend_centers[r][j].append(dends[ordered_pos[0][0]])
        dend_centers[r][j].append(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])

        
with open(pickle_file, 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump(irs, f)
f.close()
''' Pickle file data:
    
    irs : resting input resistance for each dendrite for each potassium channel block (Control, ACh, Ka, Km, Kca, Kir) 

'''
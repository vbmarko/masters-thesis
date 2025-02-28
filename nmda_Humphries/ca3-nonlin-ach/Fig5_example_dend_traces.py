# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 12:50:24 2022

@author: rh17872 - Rachel Humphries

Runs simulation for six example dendrites - takes ~ 2 hours to run.
Plots neuron morphology and example traces from Figure 5A+B of the paper
Plots nonlinear curves for the 6 example dendrites (Figure 5C).

Rachel Humphries, Jack R. Mellor, Cian O'Donnell, 
Acetylcholine Boosts Dendritic NMDA Spikes in a CA3 Pyramidal Neuron Model,
Neuroscience, 2021,ISSN 0306-4522
https://doi.org/10.1016/j.neuroscience.2021.11.014

"""
import numpy as np
from neuron import h
from neuron.units import ms, mV
h.load_file('stdrun.hoc')
import matplotlib
import matplotlib.pyplot as plt
import ca3_synapse_functions as csf
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.gridspec as gridspec
import pickle
import add_scalebar_func as scale
plt.rcParams.update({'font.size': 14})
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.rcParams["savefig.dpi"] = 600

'''Experiment parameters'''
pickle_file = "output_variables.pkl"
experiment = ["Control", "ACh", "Ka" r'$\downarrow$', "Km" r'$\downarrow$',"Kca" r'$\downarrow$',"Kir" r'$\downarrow$']
save_exps = ["Control", "ACh", "Ka", "Km", "Kca", "Kir"]
regions = ["SR", "SLM"]
example_dends = [[12,36,27],[4,3,20]] # dendrite index for SR and SLM for example dendrites

'''Synapse parameters'''
syn_ISIs = 1 
syn_density = 1 #ums between each synapse 0,1,2,5
syn_NAratios = [0.5] #AMPA to NMDA ratio
nmda_times = [5,16]
ampa_times = [0.5,1.5] #rise and decay times 0.4,4.1
total_cond_ac = 0.0035*0.5
total_cond_pp = 0.0035*0.5#0.00064

a_ac_weights, n_ac_weights = csf.calculate_synapse_weights(syn_NAratios,total_cond_ac, "AC")
a_pp_weights, n_pp_weights = csf.calculate_synapse_weights(syn_NAratios,total_cond_pp, "PP")

for i,x in enumerate(a_ac_weights) :
    print ("AC AMPA:NMDA ratios :", x/n_ac_weights[i] )
    print ("PP AMPA:NMDA ratios :", a_pp_weights[i]/n_pp_weights[i])

'''Stimulation parameters'''
step = 1
total_syns = 20
num_synapses = range(step,total_syns+1,step)
stim_num = 1 #number of stims
stim_start = 810 #start of stim (ms)
stim_noise = 0 #noise of stim
stim_interval = 20 #interval between stims (20ms = 50hz)

'''Optimiser parameters'''
def setup_neuron_model() :
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
    #h.rm_val = rm
    
    h.gkdr_val = gkdr
    h.ghd_val = gih
    
    h.km_sh_val = km_sh
    h.ka_sh_val= ka_sh
    h.kir_sh_val= kir_sh
    h.ih_sh_val = ih_sh
    
    h.dt = 1.
    dt = 1.

    return gka, gkm, gkca, gkir, ka_ach_block, km_ach_block, kca_ach_block, kir_ach_block, dt

gka, gkm, gkca, gkir, ka_ach_block, km_ach_block, kca_ach_block, kir_ach_block, dt = setup_neuron_model()

'''Setup potassium channel conductances for each run'''
def setup_potassium_channels(save_exps,pot_channel,pot_conductance,pot_block) :
    if pot_channel in save_exps :
        pot_cond_list= [pot_conductance]*(len(save_exps)) #0.02 / 0.0068 /0.8 (*40)
        pot_cond_list = [0. if i==save_exps.index(pot_channel) else x for i,x in enumerate(pot_cond_list)]
        pot_cond_list = [pot_conductance*pot_block if save_exps[i]=="ACh" else x for i,x in enumerate(pot_cond_list)]
    else :
        pot_cond_list = [pot_conductance]*(len(save_exps))
        pot_cond_list = [pot_conductance*pot_block if save_exps[i]=="ACh" else x for i,x in enumerate(pot_cond_list)]
    return pot_cond_list
    
gkas = setup_potassium_channels(save_exps,"Ka",gka,1.) #ka_ach_block has to be set during simulation in the hoc file
gkms = setup_potassium_channels(save_exps,"Km",gkm,km_ach_block)
gkcas = setup_potassium_channels(save_exps,"Kca",gkca,kca_ach_block)
gkirs = setup_potassium_channels(save_exps,"Kir",gkir,kir_ach_block)

print ("Potassium conductances: ", gkas, gkms, gkcas, gkirs)

'''Divide and setup dendrites'''
sl_dends = []
sr_dends = []
slm_dends = []
h.distance(sec=h.soma[0]) #set the origin at the soma
sl_sr_bound = 150
sr_slm_bound = 400
rel = 0.5
sl_dends, sr_dends, slm_dends = csf.divide_dends(rel,sl_sr_bound,sr_slm_bound,sl_dends,sr_dends,slm_dends)
so_dends = range(int(h.numbasal)) #basal dendrites
all_dends = [so_dends, sl_dends, sr_dends, slm_dends]

'''Identify dendrite sections to use in simulation - either example dends list or all dendrites > 20um in length'''
so_dend_200 = []
for d in so_dends :
    if h.dendrite[d].L > 20 :
        so_dend_200.append([d]) 
sr_dend_200 = [[d] for d in sr_dends if h.apical_dendrite[d].L > 20] 
slm_dend_200 = [[d] for d in slm_dends if h.apical_dendrite[d].L > 20]
    
if len(example_dends) > 0 :
    sr_dend_200 = [sr_dend_200[i] for i in example_dends[0]]
    slm_dend_200 = [slm_dend_200[i] for i in example_dends[1]]

print ("SR dends: ", sr_dend_200, "SLM dends: ", slm_dend_200)

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
        #calculate where to position synapses on dendrite section.
        if regions[r] == "SO" :
            pos_per_dend, mid_dend,mid_syn_pos,ordered_pos = csf.synapse_positioning(h.dendrite,dends,syn_density,num_synapses[-1])    
        else :
            pos_per_dend, mid_dend,mid_syn_pos,ordered_pos = csf.synapse_positioning(h.apical_dendrite,dends,syn_density,num_synapses[-1])    
        for s,nar in enumerate(experiment) : #loops through different potassium channel blocks.
            print ("Experiment: ", experiment[s])
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
            
            for n in num_synapses : #loops through number of synapses to add
                print ("Combined synapses" , n) # adds both AMPA and NMDA inputs
                glu_syns = []
                glu_netstim = []
                glu_netcon = []
                for ns in range(n) : #loops through each synapse to add
                    d = dends[ordered_pos[ns][0]]
                    p = pos_per_dend[ordered_pos[ns][0]][ordered_pos[ns][1]]
                    #print "Synapses: ", d, p
                    if regions[r] == "SO" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_nmda_baker(h.dendrite[d],p,n_ac_weights[0],nmda_times[0],nmda_times[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.dendrite[d],p,a_ac_weights[0],ampa_times[0],ampa_times[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                    if regions[r] == "SR" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_nmda_baker(h.apical_dendrite[d],p,n_ac_weights[0],nmda_times[0],nmda_times[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.apical_dendrite[d],p,a_ac_weights[0],ampa_times[0],ampa_times[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                    if regions[r] == "SLM" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_nmda_baker(h.apical_dendrite[d],p,n_pp_weights[0],nmda_times[0],nmda_times[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.apical_dendrite[d],p,a_pp_weights[0],ampa_times[0],ampa_times[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)

                soma_v = h.Vector().record(h.soma[0](0.5)._ref_v,dt) #record soma voltage
                #record stimulated dendrite voltage
                if regions[r] == "SO" :
                    dend_v = h.Vector().record(h.dendrite[dends[ordered_pos[0][0]]](pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])._ref_v,dt)
                else:
                    dend_v = h.Vector().record(h.apical_dendrite[dends[ordered_pos[0][0]]](pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])._ref_v,dt)
                    
                t = h.Vector().record(h._ref_t,dt) #record time
                
                h.fig9b() #calls the function in the hoc file to run the simulation
                
                soma_v = np.array(soma_v)
                dend_v = np.array(dend_v)
            
                soma_volts[r][j][s].append(soma_v)
                dend_volts[r][j][s].append(dend_v)
                soma_peaks[r][j][s].append(max(soma_v[int(stim_start/dt):]))
                dend_peaks[r][j][s].append(max(dend_v[int(stim_start/dt):]))
            
            freq=0
            rin= h.Impedance() # does not measure resting Rin here
            rin.compute(freq,1)
            if regions[r] == "SO" :
                input_imp = rin.input(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]],sec=h.dendrite[dends[ordered_pos[0][0]]])        
            else:
                input_imp = rin.input(pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]],sec=h.apical_dendrite[dends[ordered_pos[0][0]]])
            irs[r][j].append(input_imp)
            
            '''Repeat simulation with AMPA only synapses'''
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
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.dendrite[d],p,a_ac_weights[0],ampa_times[0],ampa_times[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                    if regions[r] == "SR" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.apical_dendrite[d],p,a_ac_weights[0],ampa_times[0],ampa_times[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)
                    if regions[r] == "SLM" :
                        glu_syns,glu_netstim,glu_netcon = csf.insert_ampar(h.apical_dendrite[d],p,a_pp_weights[0],ampa_times[0],ampa_times[1],glu_syns,glu_netstim,glu_netcon,stim_num,stim_interval,stim_start,stim_noise)

                soma_v = h.Vector().record(h.soma[0](0.5)._ref_v,dt)
                if regions[r] == "SO" :
                    dend_v = h.Vector().record(h.dendrite[dends[ordered_pos[0][0]]](pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])._ref_v,dt)
                else:
                    dend_v = h.Vector().record(h.apical_dendrite[dends[ordered_pos[0][0]]](pos_per_dend[ordered_pos[0][0]][ordered_pos[0][1]])._ref_v,dt)
                    
                t = h.Vector().record(h._ref_t,dt)

                h.fig9b()
                
                soma_v = np.array(soma_v)
                dend_v = np.array(dend_v)
            
                ampa_soma_volts[r][j][s].append(soma_v)
                ampa_dend_volts[r][j][s].append(dend_v)
                ampa_soma_peaks[r][j][s].append(max(soma_v[int(stim_start/dt):]))
                ampa_dend_peaks[r][j][s].append(max(dend_v[int(stim_start/dt):]))
        
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

'''Save data to pickle file'''
with open(pickle_file, 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([all_dends_200,soma_volts, dend_volts, dt, time, irs, soma_peaks, dend_peaks, nmda_currs, ampa_soma_peaks,ampa_dend_peaks,ampa_soma_volts, ampa_dend_volts,dend_distances,sr_dend_200,slm_dend_200,dend_diams,so_dends,sl_dends,sr_dends,slm_dends], f)
f.close()

'''Load data from pickle file''' # If have run and saved all data previously can load it from pickle file here.

# pickle_output = "ach_exp_output_v1.pkl"

# with open("ach_exp_output_v1.pkl", 'rb') as f:  # Python 3: open(..., 'rb')
#     all_dends_200,soma_volts, dend_volts, dt, time, irs, soma_peaks, dend_peaks, nmda_currs, ampa_soma_peaks,ampa_dend_peaks,ampa_soma_volts, ampa_dend_volts,dend_distances,sr_dend_200,slm_dend_200,dend_diams,so_dends,sl_dends,sr_dends,slm_dends = pickle.load(f)
# f.close()

''' 
Pickle file data - ouput from running ach_exp.py

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


'''3D neuron points for plotting morphology'''
soma_3d = [[[]for i in range(3)]]
so_dends_3d = [[[]for i in range(3)] for i in so_dends]
sl_dends_3d = [[[]for i in range(3)] for i in sl_dends]
sr_dends_3d = [[[]for i in range(3)] for i in sr_dends]
slm_dends_3d = [[[]for i in range(3)] for i in slm_dends]
diams_3d = [[],[],[],[],[]] #for soma and each region dendrite diameters, [soma, so, sl , sr, slm]

def get_3d_points(section,dend_ind,array) : #gets 3d points from NEURON
    total = section.n3d()
    #print range(total)
    for t in range(total):
        #print t
        array[dend_ind][0].append(section.x3d(t))
        array[dend_ind][1].append(section.y3d(t))
        array[dend_ind][2].append(section.z3d(t))
    return array

soma_3d = get_3d_points(h.soma[0],0,soma_3d)
diams_3d[0].append(h.soma[0].diam)
for i,d in enumerate(so_dends) :
    so_dends_3d = get_3d_points(h.dendrite[d],i,so_dends_3d)
    diams_3d[1].append(h.dendrite[d].diam)
for i,d in enumerate(sl_dends) :
    sl_dends_3d = get_3d_points(h.apical_dendrite[d],i,sl_dends_3d)
    diams_3d[2].append(h.apical_dendrite[d].diam)
for i,d in enumerate(sr_dends) :
    sr_dends_3d = get_3d_points(h.apical_dendrite[d],i,sr_dends_3d)
    diams_3d[3].append(h.apical_dendrite[d].diam)
for i,d in enumerate(slm_dends) :
    slm_dends_3d = get_3d_points(h.apical_dendrite[d],i,slm_dends_3d)
    diams_3d[4].append(h.apical_dendrite[d].diam)

flat_dends = [[item for sublist in reg for item in sublist] for reg in all_dends_200]

def which_color(dend,region) : # function for setting color based on dendrite distance from soma
    normcolor = matplotlib.colors.Normalize(vmin=150, vmax=650)
    dend_dist = (dend_distances[region][dend][0])
    linecolor = cm.plasma_r(normcolor(dend_dist))
    return linecolor, normcolor

'''Plot morphology using 3D co-ordinates''' 
fig = plt.figure(figsize=(7,6))
ax1 = plt.axes(projection='3d')
for i,s in enumerate(soma_3d) :
    ax1.plot(xs=s[0],ys=s[2],zs=s[1], color = 'lightgrey', linewidth=diams_3d[0][i])
for i,d in enumerate(so_dends_3d) :
    ax1.plot(xs=d[0],ys=d[2],zs=d[1],color="lightgrey", linewidth=diams_3d[1][i])
for i,d in enumerate(sl_dends_3d) :
    ax1.plot(xs=d[0],ys=d[2],zs=d[1],color = "lightgrey", linewidth=diams_3d[2][i])
for i,d in enumerate(sr_dends_3d) :
    if sr_dends[i] in flat_dends[0]: 
        linecolor,normcolor = which_color(flat_dends[0].index(sr_dends[i]),0) # 0 or SR
        ax1.plot(xs=d[0],ys=d[2],zs=d[1],color=linecolor, linewidth=diams_3d[3][i])
    else:
        ax1.plot(xs=d[0],ys=d[2],zs=d[1],color="lightgrey", linewidth=diams_3d[3][i])
for i,d in enumerate(slm_dends_3d) :
    if slm_dends[i] in flat_dends[1]: 
        linecolor,normcolor = which_color(flat_dends[1].index(slm_dends[i]),1) # 0 or SR
        ax1.plot(xs=d[0],ys=d[2],zs=d[1],color=linecolor , linewidth=diams_3d[4][i])
    else:
        ax1.plot(xs=d[0],ys=d[2],zs=d[1],color="lightgrey", linewidth=diams_3d[4][i])
cbar = fig.colorbar(cm.ScalarMappable(norm=normcolor,cmap='plasma_r'), ax=ax1, shrink = 0.75, aspect=15, pad=-0.25, orientation="horizontal")
cbar.set_label('Dendrite distance from soma (\u03bcm)', rotation=0, labelpad=5)
ax1.view_init(elev=0., azim=46.)
ax1.set_axis_off()
ax1.autoscale(tight=True)
ax1.set_xlim3d(-220.725,100) #-171.725,100
ax1.set_ylim3d(-220.3,100) #-159.3, 100
ax1.set_zlim3d(-180,310)
plt.savefig("ca3_recon.png", dpi=600)

'''Example dendrites''' # only need this if have run more than the 6 example dendrites
# regions = ["SR", "SLM"]
# flatten_distances = [[item for sublist in reg for item in sublist] for reg in dend_distances]
# ordered_distances = [sorted(reg) for reg in flatten_distances]
# dists = [[3,len(flatten_distances[0])//2-5,-3] ,[2,len(flatten_distances[1])//2,-2]]
# example_dends = [[] for r in regions]
# for i,r in enumerate(dists) :
#     for d in r :
#         example_dends[i].append(flatten_distances[i].index(ordered_distances[i][d]))

'''Plot example traces (from Figure 5b)'''

eg_weights = [1,9,16]
stim_start = 810 #start of stim (ms)
experiment = ["Control", "ACh", "Ka" r'$\downarrow$', "Km" r'$\downarrow$',"Kca" r'$\downarrow$',"Kir" r'$\downarrow$']
save_exps = ["Control", "ACh", "Ka", "Km", "Kca", "Kir"]

def plot_eg_voltage_traces(ax,volts,i,j,a,exp,ymax_lim,ymin_lim,scale_time,scale_mv,xmax) :
    '''
    ax : axes to plot
    volts : voltage array
    i : 0 or 1 (SR or SLM dendrite)
    j : dendrite index for accessing voltages (in this simulation example just need 'a' parameter)
    a : example dendrite index (0-2)
    exp : 0-5 index from experiment list
    ymax_lim : y-max limit
    ymin_lim : y-min limit
    scale_time : scalebar x size (ms)
    scale_mv : scalebar y size (mv)
    xmax : ms time after stim_start to plot
    '''
    for k,s in enumerate(volts[i][a][exp]):#[::2] :
        linecolor,normcolor = which_color(a,i) # 0 or SR
        if k in eg_weights :
            thick, = ax.plot(time[int((stim_start-10)/dt):int((stim_start+xmax)/dt)],s[int((stim_start-10)/dt):int((stim_start+xmax)/dt)],linewidth=1.5,color = linecolor)
        else:
            ax.plot(time[int((stim_start-10)/dt):int((stim_start+xmax)/dt)],s[int((stim_start-10)/dt):int((stim_start+xmax)/dt)],linewidth=0.5,color = linecolor)
    if i ==1 and a==2 :
        ax.annotate(experiment[exp],(0.5,1.1), xycoords='axes fraction', ha='center', va='center', size=10)
    ax.set_ylim(ymin_lim,ymax_lim)
    if i ==1 and a==2 and exp==1 :            
        scale.scale_bar(ax,scale_time,scale_mv,str(scale_time)+"ms","",1, barwidth=1, sepx=3, sepy=-1.) #ax,sizex,sizey,labelx,labely,location
    elif exp == 0 :
        ax.yaxis.set_visible(True) 
        ax.set_frame_on(True)                
        ax.spines['bottom'].set_visible(False)
        ax.xaxis.set_visible(False)
        #ax.spines['left'].set_visible(True)
    else:    
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.set_frame_on(False)

colors = ["black","royalblue","deepskyblue","darkorange","lightseagreen","orangered"]

'''Calculate NMDA-mediated voltage by subtracting AMPA-only voltage'''
nmda_volts_dend = [[[[] for i in experiment] for j in x]for x in all_dends_200]
nmda_amps_dend = [[[[] for i in experiment] for j in x]for x in all_dends_200]

for i,r in enumerate(soma_volts) :
    for j,dend in enumerate(r) : 
        for k, dens in enumerate(dend) : #loops through experiments
            for l, ns in enumerate(dens) : #loops through synapse numbers 
                nm_volts_dend = [x1 - x2 for (x1, x2) in zip(dend_volts[i][j][k][l],ampa_dend_volts[i][j][k][l])] # subtract AMPA-only voltages
                nmda_volts_dend[i][j][k].append(nm_volts_dend) 
                nmda_amps_dend[i][j][k].append(max(nm_volts_dend[int((stim_start-2)/dt):])) #amplitude of NMDA response

def plot_weight_vs_amps_fig5(ax,x,y,i,j,a,exps_to_plot,x_label) :
    '''
    Function plots no.of synapses vs NMDA voltage amplitude (Fig 5c plots)
    ax : axes to plot
    x : x-values = number of synapses
    y : y-values = NMDA response amplitudes 
    i : 0 or 1 (SR or SLM dendrite)
    j : dendrite index for accessing voltages (in this simulation example just need 'a' parameter)
    a : example dendrite index (0-2)
    exps_to_plot : 0-5 index from experiment list
    x_label : x-axis label
    '''
    for k,exp in enumerate(exps_to_plot) :
        ax.plot(x,y[i][a][exp],"-",label=experiment[exp], color = colors[exp], linewidth=1) #plot soma
    if i ==0 and a==0 :
        ax.set_xlabel(x_label)
    else:
        ax.set_xticklabels([])
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0,ymax=70)
    if i ==1 and a==2 :
        ax.legend()

save_exps = ["Control", "ACh", "Ka", "Km", "Kca", "Kir"] #experiment simulations

'''Can choose which 2 experiment simulations to plot'''

exps_to_plot = ["Control","ACh"] #choose which experiments to plot in Figure 5b+c

exp_index = [save_exps.index(a) for a in exps_to_plot]

plt.rcParams.update({'font.size': 10})
plt.rc('ytick', labelsize=8)    # fontsize of the tick labels

fig = plt.figure(figsize = (5,8))
gs1 = fig.add_gridspec(nrows=6, ncols=2,left=0.12,right = 0.60, wspace=0.0005, hspace = 0.1) 
gs2 = fig.add_gridspec(nrows=6, ncols=1, left=0.72, right=0.99, wspace=0.0, hspace = 0.2)

for i,r in enumerate(example_dends) : 
    for a,j in enumerate(r) :
        total = (len(example_dends)*len(r)) - 1
        ax = fig.add_subplot(gs1[total- ((i*len(r))+a),0]) #control
        plot_eg_voltage_traces(ax,dend_volts, i,j,a,exp_index[0], 0,-90,20,0,70)
        ax = fig.add_subplot(gs1[total- ((i*len(r))+a),1]) #ACh
        plot_eg_voltage_traces(ax,dend_volts, i,j,a,exp_index[1],0,-90,20,0,70)
        ax = fig.add_subplot(gs2[total- ((i*len(r))+a),0])        
        plot_weight_vs_amps_fig5(ax,num_synapses,nmda_amps_dend,i,j,a,exp_index,"Number of\nsynapses") 

fig.text(0.63, 0.5, 'Amplitude of NMDA voltage response (mV)', va='center', rotation='vertical', size=10)
fig.text(0.01, 0.5, 'Voltage (mV)', va='center', rotation='vertical', size=10)
plt.savefig("example_traces.png",dpi=600)
plt.show()






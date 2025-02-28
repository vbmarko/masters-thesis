# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 11:29:43 2018

@author: rh17872

CA3 single neuron model functions 

"""
from neuron import h
import numpy as np

'''Calculates synapse weights based on AMPA:NMDA ratios'''
def calculate_synapse_weights(syn_NAratios,total_cond,synapse) :
    n_weights = []
    a_weights = []
    for x in syn_NAratios :
        a_weights.append((total_cond/(x+1))*x)
        n_weights.append(total_cond - ((total_cond/(x+1))*x))
    print (synapse, "weights: ", a_weights, n_weights)
    return a_weights,n_weights

'''Divides apical dendrites into SL, SR, SLM regions based on distance from soma'''
def divide_dends(rel,sl_sr_bound,sr_slm_bound,sl_dends,sr_dends,slm_dends) :
    for dend in range(int(h.numapical)) : #loop through apical dendrites

        if h.distance(rel,sec=h.apical_dendrite[dend]) < sl_sr_bound : #if distance from soma is closer that SL/SR boundary
            sl_dends.append(dend) #add to SL dends
        if h.distance(rel,sec=h.apical_dendrite[dend]) >= sl_sr_bound and h.distance(rel,sec=h.apical_dendrite[dend]) < sr_slm_bound :
            sr_dends.append(dend)            
        if h.distance(rel,sec=h.apical_dendrite[dend]) >= sr_slm_bound : #if distance from soma is further than SLM boundary
            slm_dends.append(dend)  
    
    return sl_dends, sr_dends,slm_dends

'''Calculate dendrite diameters'''
def get_dend_diameters(all_dends_200, regions) :
    dend_diams = [[[]for j in i]for i in all_dends_200]
    for i,region in enumerate(all_dends_200) :
        for j,group in enumerate(region) :
            #print regions[i], group, print 
            for dend in group :
                if regions[i] == "SO" :
                    dend_diams[i][j].append(h.dendrite[dend].diam)
                    print ("SO dend: ", dend, h.dendrite[dend].diam)
                else :
                    print (regions[i], "dend: ",dend, h.apical_dendrite[dend].diam )
                    dend_diams[i][j].append(h.apical_dendrite[dend].diam)
    return dend_diams

'''Function for positioning synapses evenly along dendrite section - starting from the centre'''
def synapse_positioning(apical_basal,dends,density,num_synapses) :
    total_length = 0
    lengths = []
    length_pos = []
    positions_per_dend = [[] for i in dends]
    for d in dends : #loop through dendrites
        total_length+=apical_basal[d].L #calculate the total length of all sections
        lengths.append(apical_basal[d].L) #append the lengths to a list
        length_pos.append(total_length) #append the length of each dend after being added
    midpoint = total_length/2. #calculate midpoint
    startpoint = midpoint - ((density*num_synapses)/2.) #calculate startpoint
    endpoint = startpoint + (density*num_synapses) #calculate endpoint
    #print midpoint, startpoint,endpoint
    all_pos = np.linspace(startpoint,endpoint,num=num_synapses) #create a list of all positions
    syn_order = [int(num_synapses/2)] #[10]
    for i,n in enumerate(range((num_synapses)-1)) :
        if (i % 2) == 0 :
            syn_order.append(syn_order[-1] - int(i+1))
        else:
            syn_order.append(syn_order[-1] + int(i+1))            
    #print ("Syn order: ", syn_order)
    #print ("Syn pos (", str(num_synapses), " synapses) :", all_pos)
    for i,l in enumerate(length_pos): #loops through end of each dendrite
        for j,p in enumerate(all_pos) : #loops through positions of synapses
            if i == 0 :
                if p < l :
                    pos = p/lengths[i]
                    positions_per_dend[i].append(pos)
                    if j == len(all_pos)/2 :
                        mid_dend = dends[i]
                        mid_syn_pos = pos
            else :
                if p < l and p >= length_pos[i-1] :
                    #print p, length_pos[i-1] 
                    pos = (p-length_pos[i-1])/lengths[i]
                    positions_per_dend[i].append(pos)
                    #print p-length_pos[i-1]
                    #print (p-length_pos[i-1])
                    if j == len(all_pos)/2 :
                        mid_dend = dends[i]
                        mid_syn_pos = pos
    #print ("Syn pos (", str(num_synapses), " synapses) :", positions_per_dend)   
    pos_indices = []
    for i,d in enumerate(positions_per_dend) :
        for j,p in enumerate(d) :
            pos_indices.append([i,j])
    #print ("Indices: ",pos_indices)
    ordered_pos_indices = []
    for i in syn_order :
        ordered_pos_indices.append(pos_indices[i])
    #print ("Ordered indices :", ordered_pos_indices)
    return positions_per_dend, mid_dend, mid_syn_pos, ordered_pos_indices


'''Insert NMDA synapse'''
def insert_nmda_baker(apical_basal,pos,weight,tauR,tauD,syn_list,netstim_list,netcon_list,stim_no,stim_interval,stim_start,stim_noise) :
    syn = h.Exp2NMDAR(apical_basal(pos))
    syn_list.append(syn)
    syn.tau1 = tauR
    syn.tau2 = tauD
    netstim= (h.NetStim(pos))
    netstim.number= stim_no
    netstim.interval= stim_interval
    netstim.start= stim_start
    netstim.noise = stim_noise 
    netstim_list.append(netstim)
    netcon=h.NetCon(netstim,syn,0,0,weight)
    netcon_list.append(netcon)
    return syn_list,netstim_list,netcon_list

'''Insert AMPA synapse'''
def insert_ampar(apical_basal,pos,weight,tauR,tauD,syn_list,netstim_list,netcon_list,stim_no,stim_interval,stim_start,stim_noise) :
    syn = h.ExpAMPAR(apical_basal(pos))
    syn_list.append(syn)
    syn.tau1 = tauR
    syn.tau2 = tauD
    netstim= (h.NetStim(pos))
    netstim.number= stim_no
    netstim.interval= stim_interval
    netstim.start= stim_start
    netstim.noise = stim_noise 
    netstim_list.append(netstim)
    netcon=h.NetCon(netstim,syn,0,0,weight)
    netcon_list.append(netcon)
    return syn_list,netstim_list,netcon_list

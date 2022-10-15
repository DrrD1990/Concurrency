#!/usr/bin/env python
# coding: utf-8

# In[1]:


import networkx as nx
import random
import math
import numpy as np
import numpy.random as rnd
from collections import defaultdict
import matplotlib.pyplot as plt


# In[2]:


class _ListDict_(object):
    r'''
    The direct Gillespie algorithm will involve a step that randomly select an element from a set. 
    
    This class will allow me to select a random element uniformly, and then either add it or delete it from  
    a set.  
    '''
    def __init__(self):
        self.item_to_position = {}
        self.items = []
    
    def __len__(self):
        return len(self.items)
    
    def __contains__(self, item):
        return item in self.item_to_position
    
    def insert(self, item):
        if self.__contains__(item):
            self.remove(item)
            
    def update(self, item):
        r'''
        If not present, then inserts the thing
        '''
        if item in self:
            return
        self.items.append(item)
        self.item_to_position[item] = len(self.items)-1
    
    def remove(self, choice):
        position = self.item_to_position.pop(choice)
        last_item = self.items.pop()
        if position != len(self.items):
            self.items[position] = last_item
            self.item_to_position[last_item] = position
            
    def choose_random(self):
        return random.choice(self.items)
    
    def random_removal(self):
        r'''uses to choose and then remove a random node'''
        choice = self.choose_random()
        self.remove(choice)
        return choice
    
    def total_weight(self):
        return len(self)


# In[3]:


def Model_3(G, H, a, b, beta, mu, tmin, tmax):
    r'''    
    
    Performs SIS simulations for epidemics.
    
    **H** networkx Graph
        The underlying aggregated network
    **G**
        The network that has same nodes as H but no edges
    **a**
        rate from node being low activity state to high activity state
    **b**
        rate from node being high activity state to low activity state    
    **beta**
        infection rate per edge
    **mu**
        recovery rate per node 
    '''
    
    activeness = defaultdict(lambda : 'N')#'N' for low activity state
    for node in G.nodes():#randomly chooses nodes in G according to the given probability and then add it to the set of high active nodes
        if random.random() < (math.sqrt(10)-3.0)/math.sqrt(10):
            activeness[node] = 'A' #'A'for high activity state
            for nbr in H.neighbors(node):#add edges being active due to the change of states of nodes
                G.add_edge(node,nbr)
                    
    a = float(a) 
    b = float(b)
    beta = float(beta) 
    mu = float(mu)     
    
    I = [G.order()]
    S = [G.order()-I[0]]
    times = [tmin]

    t = tmin

    status = defaultdict(lambda : 'S')
    for node in G.nodes():#all nodes in G are infectious at the initial time
        status[node] = 'I'
        
    infecteds = _ListDict_()    #set of infected nodes
    IS_links = _ListDict_()     #set of edges linked between infectious and susceptible nodes
    inactivateds = _ListDict_() #set of nodes with low activity state
    activateds = _ListDict_()   #set of nodes with high activity state
    
    for node in G.nodes():
        infecteds.update(node)
        for nbr in G.neighbors(node):
            if status[nbr] == 'S':
                IS_links.update((node, nbr))
    
    for node in G.nodes():
        if activeness[node]=='N':
            inactivateds.update(node)
        else:
            activateds.update(node)

    total_recovery_rate = mu*infecteds.total_weight()#mu*I_sum
    total_transmission_rate = beta*IS_links.total_weight()#beta*IS_sum
    total_node_on_rate = a*inactivateds.total_weight()#a*inactivateds_sum
    total_node_off_rate = b*activateds.total_weight()#b*activateds_sum
    total_rate = total_recovery_rate + total_transmission_rate + total_node_on_rate + total_node_off_rate
    delay = random.expovariate(total_rate)
    t += delay

    while infecteds and t<tmax:#choose the next event based on its probability
        events = ['recovery', 'infection', 'nodeoff', 'nodeon']
        rand = random.choices(events, weights = [total_recovery_rate, total_transmission_rate, total_node_off_rate, total_node_on_rate], k=1)
        if rand == ['recovery']:
            recovering_node = infecteds.random_removal()
            status[recovering_node] = 'S'
            
            
            for nbr in G.neighbors(recovering_node):#updating IS_links
                if nbr == recovering_node:
                    continue
                elif status[nbr] == 'S':
                    IS_links.remove((recovering_node, nbr))
                else:
                    IS_links.update((nbr, recovering_node))
            times.append(t)
            S.append(S[-1]+1)#one more susceptible
            I.append(I[-1]-1)#one less infected
            
        elif rand == ['infection']:
            transmitter, recipient = IS_links.choose_random()
            status[recipient] = 'I'
            infecteds.update(recipient)
            for nbr in G.neighbors(recipient):
                if status[nbr] == 'S':
                    IS_links.update((recipient, nbr))
                elif nbr != recipient:
                    IS_links.remove((nbr, recipient))     
            times.append(t)
            S.append(S[-1]-1)
            I.append(I[-1]+1)
            
        elif rand == ['nodeoff']:#one high active node being low active
            inactivating_node = activateds.random_removal()
            activeness[inactivating_node] = 'N'
            inactivateds.update(inactivating_node)
            for nbr in H.neighbors(inactivating_node):
                if activeness[nbr] == 'N' and nbr != inactivating_node:
                    G.remove_edge(inactivating_node, nbr)
                    if status[inactivating_node] == 'I' and status[nbr] == 'S':
                        IS_links.remove((inactivating_node, nbr))
                    elif status[inactivating_node] == 'S' and status[nbr] == 'I':
                        IS_links.remove((nbr, inactivating_node))           
            times.append(t)
            S.append(S[-1])
            I.append(I[-1])
            
        else:#one low active node being high active
            activating_node = inactivateds.random_removal()
            activeness[activating_node] = 'A'
            activateds.update(activating_node)
            for nbr in H.neighbors(activating_node):
                if activeness[nbr] == 'N':
                    G.add_edge(activating_node, nbr)
                    if status[activating_node] == 'I' and status[nbr] == 'S':
                        IS_links.update((activating_node, nbr))
                    elif status[activating_node] == 'S' and status[nbr] == 'I':
                        IS_links.update((nbr, activating_node))          
            times.append(t)
            S.append(S[-1])
            I.append(I[-1])
            
        total_recovery_rate = mu*infecteds.total_weight()#updating rate of each event
        total_transmission_rate = beta*IS_links.total_weight()
        total_node_on_rate = a*inactivateds.total_weight()
        total_node_off_rate = b*activateds.total_weight()
        total_rate = total_recovery_rate + total_transmission_rate + total_node_on_rate + total_node_off_rate
        if total_rate > 0:
            delay = random.expovariate(total_rate)
        else:
            delay = float('Inf')
        t += delay
    return I[-1]


# In[4]:


tmin=0
tmax=10000

a = math.sqrt(10)-3.0 #rate from node being low active to high active
b = 3.0               #rate from node being high active to low active
beta_min = 1.4
beta_max = 10.01
beta_gap = 0.25
mu = 1.0
iteration = 1000 #number of simulations for each value of beta.
beta_list = list(np.arange(beta_min, beta_max, beta_gap))


# In[5]:


fraction_of_infecteds = []#store the average of fraction of infectious nodes for each value of beta
for beta in np.arange(beta_min, beta_max, beta_gap):
    fraction_list = []
    for i in range(iteration):
        fh_2 = open("H_er.adjlist", 'rb') #assign the underlying network H
        H = nx.read_adjlist(fh_2)
        fh_2.close()
        n = H.order()
        G = nx.Graph()
        G.add_nodes_from(H.nodes())
        
        f = Model_3(G, H, a, b, beta, mu, tmin, tmax)
        fraction = f / n
        fraction_list.append(fraction)
    fraction_mean = np.mean(fraction_list)
    fraction_of_infecteds.append(fraction_mean)

np.savetxt('M3_q0.1_f_b.txt', fraction_of_infecteds)


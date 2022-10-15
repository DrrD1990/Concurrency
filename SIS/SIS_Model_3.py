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
        choice = self.choose_random()
        self.remove(choice)
        return choice
    
    def total_weight(self):
        return len(self)


# In[3]:


def Model_3(G, H, a, b, beta, mu, tmin, tmax):
    
    activeness = defaultdict(lambda : 'N')
    for node in G.nodes():
        if random.random() < (math.sqrt(10)-3.0)/math.sqrt(10):
            activeness[node] = 'A'
            for nbr in H.neighbors(node):
                G.add_edge(node,nbr)
                    
    a = float(a) #the rate of a node from inactive to active
    b = float(b) #the rate of a node from active to inactive
    beta = float(beta) #transimition rate
    mu = float(mu)     #recovery rate
    
    I = [G.order()]
    S = [G.order()-I[0]]
    times = [tmin]

    t = tmin

    status = defaultdict(lambda : 'S')
    for node in G.nodes():
        status[node] = 'I'
        
    infecteds = _ListDict_()
    IS_links = _ListDict_()
    inactivateds = _ListDict_()
    activateds = _ListDict_()
    
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

    total_recovery_rate = mu*infecteds.total_weight()
    total_transmission_rate = beta*IS_links.total_weight()
    total_node_on_rate = a*inactivateds.total_weight()
    total_node_off_rate = b*activateds.total_weight()
    total_rate = total_recovery_rate + total_transmission_rate + total_node_on_rate + total_node_off_rate
    delay = random.expovariate(total_rate)
    t += delay

    while infecteds and t<tmax:
        events = ['recovery', 'infection', 'nodeoff', 'nodeon']
        rand = random.choices(events, weights = [total_recovery_rate, total_transmission_rate, total_node_off_rate, total_node_on_rate], k=1)
        if rand == ['recovery']:
            recovering_node = infecteds.random_removal()
            status[recovering_node] = 'S'
            
            
            for nbr in G.neighbors(recovering_node):
                if nbr == recovering_node:
                    continue
                elif status[nbr] == 'S':
                    IS_links.remove((recovering_node, nbr))
                else:
                    IS_links.update((nbr, recovering_node))
            times.append(t)
            S.append(S[-1]+1)
            I.append(I[-1]-1)
            
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
            
        elif rand == ['nodeoff']:
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
            
        else:
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
            
        total_recovery_rate = mu*infecteds.total_weight()
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
tmax=100

a = math.sqrt(10)-3.0
b = 3.0
beta_min = 1.4
beta_max = 10.01
beta_gap = 0.25
mu = 1.0
iteration = 1
beta_list = list(np.arange(beta_min, beta_max, beta_gap))


# In[5]:


fraction_of_recovered = []
for beta in np.arange(beta_min, beta_max, beta_gap):
    fraction_list = []
    for i in range(iteration):
        fh_2 = open("H_er.adjlist", 'rb') #The underlying static network doesn't change when we run the simulation, so I write it here outside the loop.
        H = nx.read_adjlist(fh_2)
        fh_2.close()
        n = H.order()
        G = nx.Graph()
        G.add_nodes_from(H.nodes())
        
        f = Model_3(G, H, a, b, beta, mu, tmin, tmax)
        fraction = f / n
        fraction_list.append(fraction)
    fraction_mean = np.mean(fraction_list)
    fraction_of_recovered.append(fraction_mean)

np.savetxt('M3_q0.1_f_b.txt', fraction_of_recovered)


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
    The Laplace Gillespie algorithm will involve a step that randomly select an element from a set. 
    
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


def Model_1prime(G, H, b, beta, mu, C1, C2, lamda1, lamda2, D1, D2, tmin, tmax):
    r'''    
    
    Performs SIS simulations for epidemics.
    
    **H** networkx Graph
        The underlying aggregated network
    **G**
        The network that has same nodes as H but no edges
    **b**
        rate from edge being active to inactive
    **beta**
        infection rate per edge
    **mu**
        recovery rate per node
    **C1,C2,lamda1,lamda2**
        coefficients of the mixture of two exponential distributions
    **D1,D2**
        coefficients of the mixture of two exponential distributions for the initial waiting time 
    '''
            
    for edge in H.edges():#randomly chooses edges in H according to the given probability and then add them to G
        if random.random()<0.1:
            G.add_edge(*edge)
            
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
        
    infecteds = _ListDict_()  #set of infected nodes
    IS_links = _ListDict_()   #set of edges linked between infectious and susceptible nodes
    P_links = _ListDict_()    #set of present edges in G 
    D_links = _ListDict_()    #set of disappeared edges


    for node in G.nodes():
        infecteds.update(node)
        for nbr in G.neighbors(node):
            if status[nbr] == 'S':
                IS_links.update((node, nbr))
                
    for edge in H.edges():
        if edge in G.edges():
            P_links.update(edge)
        else:
            D_links.update(edge)
                         
    total_edge_on_rate = 0
    rates = ["r1", "r2"]
    for edge in H.edges():#assign the edge active rate for each initial inactive edge in G 
        if edge not in G.edges():
            prob = np.random.choice(rates, p = [D1, D2])
            if prob == "r1":
                H[edge[0]][edge[1]]['name'] = lamda1
            else:
                H[edge[0]][edge[1]]['name'] = lamda2
            total_edge_on_rate = total_edge_on_rate + H[edge[0]][edge[1]]['name']
   
    

    total_recovery_rate = mu*infecteds.total_weight()#mu*I_sum
    total_transmission_rate = beta*IS_links.total_weight()#beta*IS_sum
    total_edge_off_rate = b*P_links.total_weight()#a*D_sum
    total_rate = total_recovery_rate + total_transmission_rate + total_edge_on_rate + total_edge_off_rate
    delay = random.expovariate(total_rate)
    t += delay
    
    events = ["recovery", "infection", "edgeon", "edgeoff"]
    
    while infecteds and t<tmax:#choose the next event based on its probability
        rand = np.random.choice(events, p = [total_recovery_rate/total_rate, total_transmission_rate/total_rate, total_edge_on_rate/total_rate, total_edge_off_rate/total_rate])
        if rand == "recovery":
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
        
            total_edge_on_rate = total_edge_on_rate #updating rate of each event
            total_edge_off_rate = b*P_links.total_weight()
            total_recovery_rate = mu*infecteds.total_weight()
            total_transmission_rate = beta*IS_links.total_weight()
            total_rate = total_recovery_rate + total_transmission_rate + total_edge_on_rate + total_edge_off_rate
        
        elif rand == "infection":
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
            
            total_edge_on_rate = total_edge_on_rate
            total_edge_off_rate = b*P_links.total_weight()
            total_recovery_rate = mu*infecteds.total_weight()
            total_transmission_rate = beta*IS_links.total_weight()
            total_rate = total_recovery_rate + total_transmission_rate + total_edge_on_rate + total_edge_off_rate
        
        elif rand == "edgeon":#one inactive edge being active
            u, v = D_links.random_removal()
            G.add_edge(u, v)
            P_links.update((u, v))
            if status[u] == 'I' and status[v] == 'S':
                IS_links.update((u, v))
            elif status[u] == 'S' and status[v] == 'I':
                IS_links.update((v, u))
            times.append(t)
            S.append(S[-1])
            I.append(I[-1])
            
            total_edge_on_rate = total_edge_on_rate - H[u][v]['name']
            total_edge_off_rate = b*P_links.total_weight()
            total_recovery_rate = mu*infecteds.total_weight()
            total_transmission_rate = beta*IS_links.total_weight()
            total_rate = total_recovery_rate + total_transmission_rate + total_edge_on_rate + total_edge_off_rate
        
        else:#one active edge being inactive
            j, k = P_links.random_removal()
            G.remove_edge(j, k)
            D_links.update((j, k))
            if status[j] == 'I' and status[k] == 'S':
                IS_links.remove((j, k))
            elif status[j] == 'S' and status[k] == 'I':
                IS_links.remove((k, j))
            prob = np.random.choice(rates, p = [C1, C2])#updating the rate of being active again for this edge
            if prob == "r1":
                H[j][k]['name'] = lamda1
            else:
                H[j][k]['name'] = lamda2
            times.append(t)
            S.append(S[-1])
            I.append(I[-1])
            
            total_edge_on_rate = total_edge_on_rate + H[j][k]['name']
            total_edge_off_rate = b*P_links.total_weight()
            total_recovery_rate = mu*infecteds.total_weight()
            total_transmission_rate = beta*IS_links.total_weight()
            total_rate = total_recovery_rate + total_transmission_rate + total_edge_on_rate + total_edge_off_rate
        
        if total_rate > 0:
            delay = random.expovariate(total_rate)
        else:
            delay = float('Inf')
        t += delay
    
    return I[-1]


# In[4]:


tmin=0
tmax=10000

c = 7.0-2*math.sqrt(10)#rate from node being low active to high active in Model 2
d = 9*math.sqrt(10)-27.0#rate from node being high active to low active in Model 2

A = (c**2+6*c*d+d**2)/4
C1 = c*(0.5+(c-d)/(4*np.sqrt(A)))*(1/(0.5*(3*c+d)-np.sqrt(A)))
C2 = c*(0.5-(c-d)/(4*np.sqrt(A)))*(1/(0.5*(3*c+d)+np.sqrt(A)))

lamda1 = 0.5*(3*c+d)-np.sqrt(A)
lamda2 = 0.5*(3*c+d)+np.sqrt(A)

D1 = C1*lamda2/(C1*lamda2+C2*lamda1)
D2 = C2*lamda1/(C1*lamda2+C2*lamda1)


# In[5]:


b = 18*math.sqrt(10)-54.0
beta_min = 1.75
beta_max = 2.0
beta_gap = 0.15
mu = 1.0
iteration = 1000 #number of simulations for each value of beta.
beta_list = list(np.arange(beta_min, beta_max, beta_gap))


# In[6]:


fraction_of_infecteds = []#store the average of fraction of infectious nodes for each value of beta
for beta in np.arange(beta_min, beta_max, beta_gap):
    fraction_list = []
    for i in range(iteration):
        fh_2 = open("H_er.adjlist", 'rb')#assign the underlying network H
        H = nx.read_adjlist(fh_2)
        fh_2.close()
        n = H.order()
        G = nx.Graph()
        G.add_nodes_from(H.nodes())
    
        f = Model_1prime(G, H, b, beta, mu, C1, C2, lamda1, lamda2, D1, D2, tmin, tmax)
        fraction = f / n
        fraction_list.append(fraction)
    fraction_mean = np.mean(fraction_list)
    fraction_of_infecteds.append(fraction_mean)

np.savetxt('M1p_q0.1_f_a.txt', fraction_of_infecteds)


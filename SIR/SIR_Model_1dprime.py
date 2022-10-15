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


def Model_1dprime(G, H, a, beta, mu, C1, C2, lamda1, lamda2, D1, D2, initial_infecteds, initial_recovereds, tmin = 0, tmax = float('Inf')):
    r'''    
    
    Performs SIR simulations for epidemics.
    
    **H** networkx Graph
        The underlying aggregated network
    **G**
        The network that has same nodes as H but no edges
    **a**
        rate from edge being inactive to active
    **beta**
        infection rate per edge
    **mu**
        recovery rate per node
    **C1,C2,lamda1,lamda2**
        coefficients of the mixture of two exponential distributions
    **D1,D2**
        coefficients of the mixture of two exponential distributions for the initial waiting time 
    **initial_infecteds**   
        one randomly choosed node in G is initially infected
    **initial_recovereds**
        empty
    '''
    for edge in H.edges():#randomly chooses edges in H according to the given probability and then add them to G
        if random.random()<0.1:
            G.add_edge(*edge)
            
    I = [len(initial_infecteds)]
    R = [len(initial_recovereds)]
    S = [G.order()-I[0]-R[0]]
    times = [tmin]

    t = tmin

    status = defaultdict(lambda : 'S')
    for node in initial_infecteds:
        status[node] = 'I'
    for node in initial_recovereds:
        status[node] = 'R'
        
    infecteds = _ListDict_() #set of infected nodes
    IS_links = _ListDict_()  #set of edges linked between infectious and susceptible nodes
    P_links = _ListDict_()   #set of present edges in G
    D_links = _ListDict_()   #set of disappeared edges

    for node in initial_infecteds: #add the initial infectious node to the set of infected nodes
        infecteds.update(node)
        for nbr in G.neighbors(node):#add the initial IS_links to the set of IS_links
            if status[nbr] == 'S':
                IS_links.update((node, nbr))
    
    for edge in H.edges():
        if edge in G.edges():
            P_links.update(edge)
        else:
            D_links.update(edge)
            
    total_edge_off_rate = 0
    rates = ["r1", "r2"]
    for edge in G.edges(): #assign the edge inactive rate for each initial active edge in G 
        prob = np.random.choice(rates, p = [D1, D2])
        if prob == "r1":
            H[edge[0]][edge[1]]['name'] = lamda1
        else:
            H[edge[0]][edge[1]]['name'] = lamda2
        total_edge_off_rate = total_edge_off_rate + H[edge[0]][edge[1]]['name']
        

    total_recovery_rate = mu*infecteds.total_weight()#mu*I_sum
    total_transmission_rate = beta*IS_links.total_weight()#beta*IS_sum
    total_edge_on_rate = a*D_links.total_weight()#a*D_sum
    total_rate = total_recovery_rate + total_transmission_rate + total_edge_on_rate + total_edge_off_rate
    delay = random.expovariate(total_rate)
    t += delay
    
    events = ["recovery", "infection", "edgeon", "edgeoff"]
    
    while infecteds and t<tmax:#choose the next event based on its probability
        rand = np.random.choice(events, p = np.around([total_recovery_rate/total_rate, total_transmission_rate/total_rate, total_edge_on_rate/total_rate, total_edge_off_rate/total_rate], decimals = 10))
        if rand == "recovery":
            recovering_node = infecteds.random_removal()
            status[recovering_node] = 'R'
            
            for nbr in G.neighbors(recovering_node):
                if status[nbr] == 'S':
                    IS_links.remove((recovering_node, nbr))#updating IS_links
            times.append(t)
            S.append(S[-1])
            I.append(I[-1]-1)#one less infected
            R.append(R[-1]+1)#one more recovered
            
            total_edge_off_rate = total_edge_off_rate #updating rate of each event
            total_edge_on_rate = a*D_links.total_weight()
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
                elif status[nbr] == 'I' and nbr != recipient:
                    IS_links.remove((nbr, recipient))
            times.append(t)
            S.append(S[-1]-1)#one less susceptible
            I.append(I[-1]+1)#one more infected
            R.append(R[-1])
            
            total_edge_off_rate = total_edge_off_rate
            total_edge_on_rate = a*D_links.total_weight()
            total_recovery_rate = mu*infecteds.total_weight()
            total_transmission_rate = beta*IS_links.total_weight()
            total_rate = total_recovery_rate + total_transmission_rate + total_edge_on_rate + total_edge_off_rate
        
        elif rand == "edgeon":#one inactive edge being active
            u, v = D_links.choose_random()
            G.add_edge(u, v)
            D_links.remove((u, v))
            P_links.update((u, v))
            if status[u] == 'I' and status[v] == 'S':
                IS_links.update((u, v))
            elif status[u] == 'S' and status[v] == 'I':
                IS_links.update((v, u))
            prob = np.random.choice(rates, p = [C1, C2])#updating the rate of being inactive again for this edge
            if prob == "r1":
                H[u][v]['name'] = lamda1
            else:
                H[u][v]['name'] = lamda2    
            times.append(t)
            S.append(S[-1])
            I.append(I[-1])
            R.append(R[-1])
            
            total_edge_off_rate = total_edge_off_rate + H[u][v]['name']
            total_edge_on_rate = a*D_links.total_weight()
            total_recovery_rate = mu*infecteds.total_weight()
            total_transmission_rate = beta*IS_links.total_weight()
            total_rate = total_recovery_rate + total_transmission_rate + total_edge_on_rate + total_edge_off_rate
        
        else:#one active edge being inactive
            j, k = P_links.choose_random()
            G.remove_edge(j, k)
            P_links.remove((j, k))
            D_links.update((j, k))
            if status[j] == 'I' and status[k] == 'S':
                IS_links.remove((j, k))
            elif status[j] == 'S' and status[k] == 'I':
                IS_links.remove((k, j))
            times.append(t)
            S.append(S[-1])
            I.append(I[-1])
            R.append(R[-1])
            
            total_edge_off_rate = total_edge_off_rate - H[j][k]['name']
            total_edge_on_rate = a*D_links.total_weight()
            total_recovery_rate = mu*infecteds.total_weight()
            total_transmission_rate = beta*IS_links.total_weight()
            total_rate = total_recovery_rate + total_transmission_rate + total_edge_on_rate + total_edge_off_rate
        
        if total_rate > 0:
            delay = random.expovariate(total_rate)
        else:
            delay = float('Inf')
        t += delay
     
    return R[-1]


# In[4]:


c = math.sqrt(10)-3.0 #rate from node being low active to high active in Model 3
d = 3.0               #rate from node being high active to low active in Model 3
A = (d**2+6*c*d+c**2)/4
C1 = d*(0.5+(d-c)/(4*np.sqrt(A)))*(1/(0.5*(3*d+c)-np.sqrt(A)))
C2 = d*(0.5-(d-c)/(4*np.sqrt(A)))*(1/(0.5*(3*d+c)+np.sqrt(A)))

lamda1 = 0.5*(3*d+c)-np.sqrt(A)
lamda2 = 0.5*(3*d+c)+np.sqrt(A)

D1 = C1*lamda2/(C1*lamda2+C2*lamda1)
D2 = C2*lamda1/(C1*lamda2+C2*lamda1)


# In[5]:


a = 2*math.sqrt(10)-6.0
beta_min = 0.0
beta_max = 10.01
beta_gap = 0.25
mu = 1.0
iteration = 1000 #number of simulations for each value of beta.
beta_list = list(np.arange(beta_min, beta_max, beta_gap))


# In[6]:


fraction_of_recovered = []#store the average of fraction of final recovered nodes for each value of beta
for beta in np.arange(beta_min, beta_max, beta_gap):
    fraction_list = []
    for i in range(iteration):
        fh_2 = open("H_er.adjlist", 'rb')#assign the underlying network H
        H = nx.read_adjlist(fh_2)
        fh_2.close()
        n = H.order()
        G = nx.Graph()
        G.add_nodes_from(H.nodes())
        
        initial_infecteds=random.sample(G.nodes(), 1) #randomly choose the initial infected node
        initial_recovereds = []
        
        f = Model_1dprime(G, H, a, beta, mu, C1, C2, lamda1, lamda2, D1, D2, initial_infecteds, initial_recovereds, tmin = 0, tmax = float('Inf'))
        fraction = f / n
        fraction_list.append(fraction)
    fraction_mean = np.mean(fraction_list)
    fraction_of_recovered.append(fraction_mean)

np.savetxt('M1dp_q0.1_f_1k1.txt', fraction_of_recovered)


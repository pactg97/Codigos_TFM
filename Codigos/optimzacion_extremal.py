# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:47:19 2020

@author: 34625
"""

## Optimización extremal bi particion

import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import random 

def biparticion(nodos):
    n=len(nodos)
    orden=random.sample(nodos,n)
    tamaño=n//2
    conj1=[orden[i] for i in range(tamaño)]
    conj2=[orden[i] for i in range(tamaño,n)]
    return conj1,conj2


def calculo_B(A,nodos,n,m,k):
    B=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            B[i,j]=A[i,j]-((k[i]*k[j])/(2*m))
    return B

def modularity(B,c,n,m):
    mod=0
    for i in range(n):
        for j in range(n):
            if(c[i]==c[j]):
                mod+=B[i,j]
    mod=mod/(2*m)
    return mod

def mejor_movimiento(B,c,n,m):
    incremento=np.zeros(n)
    c0=[i for i in c]
    for i in range(n):
        incremento[i]=sum([B[i,j] for j in range(n) if c0[j]!=c0[i]]) - sum([B[i,j] for j in range(n) if c0[j]==c0[i]])
    max_incremento=max(incremento)
    movimientos=[i for i in range(n) if incremento[i]==max_incremento] 
    if(max_incremento/(2*m)<=0.001):
        return c0
    else:
        for i in movimientos:
            if(c0[i]==0): c0[i]=1
            else: c0[i]=0
        return c0

def opt_ext_bi(G):
    if(nx.is_connected(G)):
        A=nx.adjacency_matrix(G)
        nodos=list(G.nodes())
        n=G.order()
        m=G.number_of_edges()
        k=[G.degree(nodos[i]) for i in range(n)]
        B=calculo_B(A,nodos,n,m,k)
        nodos1,nodos2=biparticion(nodos)
        c=[0 for i in range(n)]
        for i in range(n):
            if(nodos[i] in nodos1):c[i]=0
            else: c[i]=1
        for i in range(10):
            mod0=None
            while(c!=mejor_movimiento(B,c,n,m)):
                c=mejor_movimiento(B,c,n,m)
                mod=modularity(B,c,n,m)
                if(mod0 is None or mod0 < mod):
                    mod0=mod
                    c_opt=c 
    comunidades=[]
    for i in set(c_opt):
        conj=set()
        for j in range(n):
            if(c_opt[j]==i): conj.add(j)
        comunidades.append(conj)
    return comunidades,mod

def opt_ext(G,k):
    c_inicial,_=opt_ext_bi(G)
    n=G.order()
    nodos=list(G.nodes())
    conj_comp=list(set(c_inicial))
    while(len(conj_comp)<k):
        for i in conj_comp:
            Gi=nx.Graph(list(G.edges()))
            cnoti=[j for j in range(n) if c_inicial[j]!=i]
            ci=[j for j in range(n) if c_inicial[j]==i]
            Gi.remove_nodes_from([nodos[j] for j in cnoti])
            ci2,_=opt_ext_bi(Gi)
            for j in range(len(ci)):
                c_inicial[ci[j]]=ci2[j]
        conj_comp=list(set(c_inicial))
    return c_inicial
    
G2=nx.Graph()
G2.add_nodes_from(range(34))
G2.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])
plt.figure(1)
nx.draw_shell(G2,with_labels=True)
print(opt_ext_bi(G2))
    
    
    
    
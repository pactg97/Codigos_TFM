# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 15:39:13 2020

@author: 34625
"""

import networkx as nx 
import matplotlib.pyplot as plt
import random 
import math

def preproceso(G):
    n=G.order()
    mapping=dict()
    nodos=list(G.nodes())
    for i in range(n):
        mapping.update({nodos[i]:i})
    H= nx.relabel_nodes(G, mapping)
    return H

## Algoritmo de Karger_stein es una variante del algoritmo de Karger que nos permite encontrar la bipartición
## que miniza la función corte con una mayor probabilidad que la del algoritmo original de  Karger

def partition(collection):
    if len(collection) == 1:
        yield [ collection ]
        return

    first = collection[0]
    for smaller in partition(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
        # put `first` in its own subset 
        yield [ [ first ] ] + smaller

def bipartition(n):
    return [i for i in list(partition(list(range(n)))) if len(i)==2]

def random_edge(G):
    edges=list(G.edges())
    edge=random.choice(edges)
    return edge[0],edge[1]

def contract(G,u,v):
    G=nx.contracted_nodes(G,u,v,self_loops=False)
    if(type(u) is tuple):tuple_u=u
    else: tuple_u=tuple([u])
    if(type(v) is tuple):tuple_v=v
    else: tuple_v=tuple([v])
    tupla=tuple_u+tuple_v
    G=nx.relabel_nodes(G,{u:tupla})
    return G

def Karger(G,t):
    while(G.order()>t):
        u,v=random_edge(G)
        G=contract(G,u,v)
    return G

def mincut(G):
    cut=G.number_of_edges()+1
    n=G.order()
    biparticiones=bipartition(n)
    nodos=list(G.nodes())
    for bipart in biparticiones:  
        nodosA=[nodos[i] for i in bipart[0]]
        nodosB=[nodos[i] for i in bipart[1]]
        edgesAB=[i for i in list(G.edges()) if((i[0] in nodosA and i[1] in nodosB) or (i[1] in nodosA and i[0] in nodosB))]
        cut1=len(edgesAB)
        if(cut1<cut): 
            cut=cut1
            A=()
            for i in nodosA:
                if(type(i) is tuple):A=A+i
                else: A=A+tuple([i])
            B=()
            for i in nodosB:
                if(type(i) is tuple):B=B+i
                else: B=B+tuple([i])
            biparticion=[set(A),set(B)]
    return biparticion,cut
    

def Karger_Stein(G):
    n=G.order()
    if(n<=6):
        return mincut(G)
    else:
        t=int(1+n/math.sqrt(2))+1
        G1=Karger(G,t)
        G2=Karger(G,t)
        bipart1,corte1=Karger_Stein(G1)
        bipart2,corte2=Karger_Stein(G2)
        if(corte1<corte2):return bipart1,corte1
        else: return bipart2, corte2
        
def MonteCarlo_Karger_Stein(G,K):
    cut=G.number_of_edges()+1
    for i in range(K):
        G1,cut1=Karger_Stein(G)
        if(cut1<cut): 
            cut=cut1
            G2=G1
    return G2,cut

## Primer ejemplo: Cuatro K5 unidos dos a dos
G=nx.MultiGraph()
for i in range(20):
    G.add_node(i)
for i in range(5):
    for j in range(5):
        if(i!=j):
            G.add_edge(i,j)
for i in range(5,10):
    for j in range(5,10):
        if(i!=j):
            G.add_edge(i,j)
for i in range(10,15):
    for j in range(10,15):
        if(i!=j):
            G.add_edge(i,j)
for i in range(15,20):
    for j in range(15,20):
        if(i!=j):
            G.add_edge(i,j)
G.add_edge(0,5)
G.add_edge(5,10)
G.add_edge(10,15)
G.add_edge(15,0)

## Represento gráficamente el grafo original
plt.figure(1)
nx.draw_shell(G,with_labels=True)
print("Ejemplo1")
## Imprimo por pantalla la bipartición y la función corte
print(MonteCarlo_Karger_Stein(G,1))


## Segundo Ejemplo: obtenido del artículo " M.E.J. Newman, Modularity and community structure 
## in networks, Proc. Natl. Acad. Sci. 103 (23) (2006) 8577–8582"
G1=nx.MultiGraph()
G1.add_nodes_from(range(34))
G1.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])
    
G2=nx.MultiGraph()
G2.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])
H=preproceso(G2)
    
## Represento gráficamente el grafo original
plt.figure(2)
nx.draw_shell(G1,with_labels=True)
print("Ejemplo2")
## Imprimo por pantalla la bipartición y la función corte
print(MonteCarlo_Karger_Stein(G1,1))
print(MonteCarlo_Karger_Stein(H,1))
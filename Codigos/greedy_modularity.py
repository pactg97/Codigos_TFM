# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 09:19:52 2020

@author: 34625
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from itertools import chain, combinations

def preproceso(G):
    n=G.order()
    mapping=dict()
    nodos=list(G.nodes())
    for i in range(n):
        mapping.update({nodos[i]:i})
    H= nx.relabel_nodes(G, mapping)
    return H

## Este algoritmo parte una partición total, es decir, cada vértice genera una comunidad formada unicamente
## por el mismo. Va uniendo parejas de comunidades de forma que el incremento de modularidad sea el máximo
## posible. Continuo con este proceso hasta alcanzar el número de comunidades en las que quiera particionar
## mi grafo

## Calcula la matriz B de la modularidad
def calculo_B(A,n,m,k):
    B=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            B[i,j]=A[i,j]-((k[i]*k[j])/(2*m))
    return B

## Calcula la modularidad para una matriz B y una partición c
def modularity(G,particion):
    m=G.number_of_edges()
    aristas=list(G.edges())
    mod=0
    grados={}
    for (i,j) in list(G.degree()):
        grados[i]=j
    for comunidad in particion:
        for i in comunidad:
            for j in comunidad:
                mod=mod+aristas.count((i,j))+aristas.count((j,i))-grados[i]*grados[j]/(2*m)
    mod=mod/(2*m)   
    return mod

def min_mod(G):
    m=G.number_of_edges()
    grados={}
    nodos=list(G.nodes())
    for (i,j) in list(G.degree()):
        grados[i]=j
    min_mod=sum([sum([-grados[i]*grados[j]/(2*m) for j in nodos]) for i in nodos])/(2*m)
    return min_mod

## Busco las uniones de comunidades de la partición c que generen el incremento de modularidad máximo
## En principio cada comunidad es representada por el índice mínimo que pertence a la comunidad.
## He considerado que si dos uniones generan el mismo incremento máximo  de modularidad entonces 
## se deben hacer simultáneamente, y si estas uniones tienen una comunidad común, entonces se deben
## unir las comunidades en el mismo paso.

## Voy uniendo comunidades partiendo de una partición total hasta llegar a un número de comunidades C o menor
def greedy_modularity(G):
    particion_inicial={frozenset([i]) for i in list(G.nodes)}
    particiones_mod=[]
    while(len(particion_inicial)>2):
        mod=min_mod(G)
        for i in particion_inicial:
            for j in particion_inicial:
                if(i!=j):
                    part_prueba=particion_inicial-{i}-{j}
                    part_prueba=part_prueba.union({i.union(j)})
                    if(mod<= modularity(G,part_prueba)):
                        part_opt=particion_inicial-{i}-{j}
                        part_opt=part_opt.union({i.union(j)})
                        mod=modularity(G,part_opt)
        particion_inicial=part_opt
        particiones_mod=particiones_mod+[(part_opt,mod)]
    return particiones_mod


## Segundo Ejemplo: obtenido del artículo " M.E.J. Newman, Modularity and community structure 
## in networks, Proc. Natl. Acad. Sci. 103 (23) (2006) 8577–8582"
G1=nx.Graph()
G1.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])

G2=nx.Graph()
G2.add_nodes_from(range(34))
G2.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])

## Represento gráficamente el grafo original
plt.figure(2)
nx.draw_shell(G2,with_labels=True)
print("Ejemplo2")
## Imprimo por pantalla las comunidades y la modularidad de la partidición
plt.figure(3)
G6=nx.karate_club_graph()
mapping={}
N=G6.order()
for i in range(N):
    mapping[i]=i+1
G7=nx.relabel_nodes(G6,mapping)
particiones=greedy_modularity(G6)
mod_max=max([j for (i,j) in particiones])
part_max=[i for (i,j) in particiones if(j==mod_max)][0]
print(part_max)
colores=["red","green","lightgreen","lightcoral","brown","white"]
colors=["red" for i in range(N)]
for k in range(len(list(part_max))):
    for i in list(part_max)[k]:
        colors[i]=colores[k]
plt.figure(figsize=(8, 6))
nx.draw_shell(G7,with_labels=True,node_color=colors,node_size=700,font_size=18)
plt.savefig('zachary_greedy.jpg')


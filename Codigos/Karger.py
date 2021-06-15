# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 09:53:25 2020

@author: 34625
"""

## Algoritmo de Karger para buscar el corte mínimo (Bipartición que minimice la función de corte)
## Consiste en contraer aristas de forma aleatoria hasta obtener unicamente dos nodos en nuestro grafo,
## de forma, que estos dos nodos están representados por nodos originales que van formar cada comunidad.

## Si nuestro grafo tiene una estructura de comunidades, obtendremos las comunidades a través del algoritmo
## con una mayor probabilidad. Es por esto que realizaré varias veces el algoritmo (Simulación de MonteCarlo)
## quedandome con la bipartición que minimice la función corte.

## Como veremos en los ejemplos si existen nodos con grado 1, el valor de corte mínimo se alcanza al aislar
## este nodo. Debido a esto en muchos casos es mejor minimizar la función de corte normalizada.

import networkx as nx 
import matplotlib.pyplot as plt
import random 

def preproceso(G):
    n=G.order()
    mapping=dict()
    nodos=list(G.nodes())
    for i in range(n):
        mapping.update({nodos[i]:i})
    H= nx.relabel_nodes(G, mapping)
    return H

## Selección de arista a contraer
def random_edge(G):
    edges=list(G.edges())
    edge=random.choice(edges)
    return edge[0],edge[1]

## Se contrae la arista uv eliminando blucles de un nodo consigo mismo
def contract(G,u,v):
    G=nx.contracted_nodes(G,u,v,self_loops=False)
    if(type(u) is tuple):tuple_u=u
    else: tuple_u=tuple([u])
    if(type(v) is tuple):tuple_v=v
    else: tuple_v=tuple([v])
    tupla=tuple_u+tuple_v
    G=nx.relabel_nodes(G,{u:tupla})
    return G

## Función corte para una biparticion A,B
def Cut(G,A,B):
    cut=0
    for i in A:
        for j in B:
            if(i in list(G.neighbors(j))): cut+=1
    return cut

## Algoritmo de Karger
def Karger(G):
    while(G.order()>2):
        u,v=random_edge(G)
        G=contract(G,u,v)
    return G,G.number_of_edges()

## Simulo K-veces el algoritmo Karger y me quedo con el corte mínimo
def MonteCarlo_Karger(G,K):
    cut=G.number_of_edges()
    for i in range(K):
        G1,cut1=Karger(G)
        if(cut1<cut): 
            cut=cut1
            G2=G1
    nodos=list(G2.nodes())
    for i in range(2):
        if(type(nodos[i]) is not tuple): nodos[i]=tuple([nodos[i]])
    return [set(i) for i in nodos],cut

## En estos casos es necesario definir los grafos como Multigrafos para que existan aristas repetidas, de esta 
## forma habrá una mayor probabilidad de contraer aristas que se encuentren dentro de una comunidad.

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
## Imprimo por pantalla la bipartición y la función corte
print("Ejemplo1")
print(MonteCarlo_Karger(G,5))

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
print(MonteCarlo_Karger(G1,10))
print(MonteCarlo_Karger(H,10))

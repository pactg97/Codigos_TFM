# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 09:21:06 2020

@author: 34625
"""
## Algoritmo Girvan-Newman hasta obtener dos componentes conexas. Se podrían cambiar los criterios
## para que siga quitando aristas hasta un número de componentes máxima o hasta quitar un número de aristas 
## máxima. También sería interesante si el conjunto de aristas a quitar ("edges_to_remove") es muy grande
## que no las quitara y para el algoritmo en ese paso antes de quitarlas, ya que puede dañar la estructura
## en comunidades. A parte de mostrar las comunidades finales, calculo la modularidad para estas particiones y
## la muestro

## Librerías utilizadas
import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np

def preproceso(G):
    n=G.order()
    mapping=dict()
    nodos=list(G.nodes())
    for i in range(n):
        mapping.update({nodos[i]:i})
    H= nx.relabel_nodes(G, mapping)
    return H

## Calculo la matriz B utilizada para el cálculo de la modularidad donde:
## B_{ij}= A_{ij}-k(i)k(j)/2m, donde k(i) es el grado del nodo i, A la matriz de adyacencia y m el 
## número de aristas
def calculo_B(A,n,m,k):
    B=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            B[i,j]=A[i,j]-((k[i]*k[j])/(2*m))
    return B

## Calculo la modularidad del grafo para una partición dada por c
def modularity(B,c,n,m):
    mod=0
    for i in c:
        for j in i:
            for k in i:
                mod+=B[j,k]
    mod=mod/(2*m)
    return mod
  
## Cálculo el valor máximo de los "betweeness" de cada arista del grafo
## También obtengo las aristas que tienen dicho valor máximo
def bet_max(g):
    d1 = nx.edge_betweenness_centrality(g) 
    list_of_tuples = list(d1.items())  
    lista_ordenada=sorted(list_of_tuples, key = lambda x:x[1],reverse=True)
    _,valor_max=lista_ordenada[0]
    # Will return in the form (a,b) 
    return lista_ordenada,valor_max
  
## Elimino de mi grafo las aristas con un valor máximo de "betweeness"
def edges_to_remove(g):
    lista,valor_max=bet_max(g)
    edges=[]
    for i in lista:
          edge,bet=i
          if(bet==valor_max): 
              edges.append(edge)
          else:
              break
    return edges
    
## Voy eliminando las aristas con "betweeness" máximo hasta que obtenga un grafo no conexo donde sus 
## componentes serán las comunidades. Cada vez que elimino las aristas con betweeness maximo debo volver 
## a calcular el valor de este para cada una de las aristas en nuestro nuevo gtafo.
    
def girvan(g): 
    a = nx.connected_components(g) 
    lena = len(list(a)) 
    while (lena <=1): 
        # We need (a,b) instead of ((a,b)) 
        edges= edges_to_remove(g) 
        for i in edges:
            u,v=i
            g.remove_edge(u,v)    
        a = nx.connected_components(g) 
        lena=len(list(a)) 
        comunidades=sorted(nx.connected_components(g), key = len, reverse=True)
    return comunidades
  
## Vamos a ver los resultados para unos ejemplos sencillos

## Primer ejemplo: Dos K4 unidos por una arista
G1=nx.Graph([(0,1),(0,2),(1,2),(2,3),(3,4),(4,5),(3,5),(6,5),(6,4),(6,3),(7,2),(7,0),(7,1)]) 

plt.figure(1)
## Represento el grafo original
nx.draw_shell(G1,with_labels=True)

## Aparece por pantalla las comunidades detectadas y la modularidad de la partición
print("Ejemplo1")
print(girvan(G1))

plt.figure(2)
## Represento el grafo tras quitar todas las aristas
nx.draw_shell(G1,with_labels=True)

## Segundo ejemplo: Cuatro K5 unidos dos a dos
G2=nx.Graph()
for i in range(20):
    G2.add_node(i)
for i in range(5):
    for j in range(5):
        if(i!=j):
            G2.add_edge(i,j)
for i in range(5,10):
    for j in range(5,10):
        if(i!=j):
            G2.add_edge(i,j)
for i in range(10,15):
    for j in range(10,15):
        if(i!=j):
            G2.add_edge(i,j)
for i in range(15,20):
    for j in range(15,20):
        if(i!=j):
            G2.add_edge(i,j)
G2.add_edge(0,5)
G2.add_edge(5,10)
G2.add_edge(10,15)
G2.add_edge(15,0)

plt.figure(3)
## Represento el grafo original
nx.draw_shell(G2,with_labels=True)

## Aparece por pantalla las comunidades detectadas y la modularidad de la partición
print("Ejemplo2")
print(girvan(G2))

plt.figure(4)
## Represento el grafo tras quitar todas las aristas
nx.draw_shell(G2,with_labels=True)

## Tercer Ejemplo: obtenido del artículo " M.E.J. Newman, Modularity and community structure 
## in networks, Proc. Natl. Acad. Sci. 103 (23) (2006) 8577–8582"
G4=nx.Graph()
G4.add_nodes_from(range(34))
G4.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])
    
G5=nx.karate_club_graph()


    
## Represento el grafo original
plt.figure(7)
nx.draw_shell(G4,with_labels=True)
## Aparece por pantalla las comunidades detectadas y la modularidad de la partición
print("Ejemplo3")
print(girvan(G4))
## Represento el grafo tras quitar todas las aristas
plt.figure(8)
nx.draw_shell(G4,with_labels=True)
plt.figure(9)
comunidades=girvan(G5)
G6=nx.karate_club_graph()
mapping={}
N=G6.order()
for i in range(N):
    mapping[i]=i+1
G7=nx.relabel_nodes(G6,mapping)
colores=[]
for i in range(0,N):
    if(i in comunidades[0]):
        colores.append("red")
    else:
        colores.append("green")
plt.figure(figsize=(8, 6))
nx.draw_shell(G7,with_labels=True,node_color=colores,node_size=700,font_size=18)

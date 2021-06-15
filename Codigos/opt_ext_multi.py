# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 10:52:28 2020

@author: 34625
"""

import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import random 

def preproceso(G):
    n=G.order()
    mapping=dict()
    nodos=list(G.nodes())
    for i in range(n):
        mapping.update({nodos[i]:i})
    H= nx.relabel_nodes(G, mapping)
    return H

## La optimización extremal trata de maximizar la función de modularidad. Para ello parte de una partición 
## aleatoria y considera como posibles movimientos cambiar un nodo de una comunidad a otra o crear una nueva 
## comunidad formada unicamente por ese nodo. 
## Buscaremos el movimiento que incremente en mayor medida la modularidad. Continuamos haciendo esto hasta que 
## se llegue a un punto que no pueda aumentar la modularidad para ningún movimiento (óptimo local)

## Cálculo de la matriz de modularidad B
def calculo_B(A,n,m,k):
    B=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            B[i,j]=A[i,j]-((k[i]*k[j])/(2*m))
    return B

## Cálculo de la modularidad para la matriz B y la partición c
def modularity(B,c,n,m):
    mod=0
    for i in c:
        for j in i:
            for k in i:
                mod+=B[j,k]
    mod=mod/(2*m)
    return mod

## Construimos una partición aleatoria de G. Desordenamos los nodos y vamos tomando elementos por la izquierda
## a través de otros números aleatorios que indican el número de elementos de cada comunidad de la partición
def particion(G):
    n=G.order()
    orden=random.sample(range(n),n)
    conj=[]
    while(len(orden)>0):
        r=random.randint(1,len(orden))
        conj1={orden.pop(0)}
        for i in range(r-1):
            conj1.add(orden.pop(0))
        conj.append(conj1)
    return conj

## Buscamos el movimiento que maximiza el incremento de modularidad
def mejor_movimiento(B,c,n,m):
    num_com=len(c)
    mod_local=np.zeros(n)
    incremento=np.zeros([n,num_com+1])
    for i in range(num_com):
        for j in c[i]:
            mod_local[j]=sum([B[j,k] for k in c[i]])   
    for i in range(n):
        for j in range(num_com):
            incremento[i,j]=sum([B[i,k] for k in c[j]])-mod_local[i]
        incremento[i,num_com]=-mod_local[i]
    max_incremento=max([max([incremento[i,j] for i in range(n)]) for j in range(num_com+1)])
    if(max_incremento==0): return 0,[[]]
    movimientos=[[(i,j) for i in range(n) if(incremento[i,j]==max_incremento)] for j in range(num_com)]
    return max_incremento,movimientos

## Por estética transformo una lista de listas en una lista con los movimientos a realizar
def trans_lista(movimientos):
    lista=[]
    for i in movimientos:
        for j in i:
            lista.append(j)
    return lista

## Aplicamos el movimiento que incrementa en mayor medida la modularidad sobre la partición c
def mover(c,movimientos):
    c1=[set(list(l)[:]) for l in c]
    num_com=len(c1)
    for (i,j) in movimientos:
        for k in range(num_com):
            if(i in c1[k]): c1[k].remove(i)
        if(j==num_com): c1.append({i})
        else: c1[j].add(i)
    return c1

## Eliminamos el conjunto vacío que nos puede quedar con alguno de los movimientos
def depurar(c):
    num_com=len(c)
    k=-1
    for i in range(num_com):
        if(c[i]==set()): k=i
    if(k>=0): c.pop(k)
    return c

## Debido a nuestro algoritmo cálcula el óptimo local para una partición inicial, repito varias 
## veces el algoritmo y me quedo con la partición que genere la mayor modularidad
def opt_ext(G,K):
    A=nx.adjacency_matrix(G)
    n=G.order()
    m=G.number_of_edges()
    k=[G.degree(i) for i in range(n)]
    mod0=None
    for i in range(K):
        c=particion(G)
        B=calculo_B(A,n,m,k)
        maxi=1
        c0=None
        while(c0 is None or modularity(B,c0,n,m)>modularity(B,c,n,m)):
            if(c0 is not None): c=[i for i in c0]
            maxi,movimientos=mejor_movimiento(B,c,n,m)
            movimientos=trans_lista(movimientos)
            c0=mover(c,movimientos)
            c0=depurar(c0)
        mod=modularity(B,c,n,m)
        if(mod0 is None or mod>mod0): 
            mod0=mod
            c_optima=c0
    return mod0,c_optima
    
## Primer ejemplo: Cuatro K5 unidos dos a dos
G=nx.Graph()
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
## Imprimo por pantalla la bipartición y la modularidad.
print(opt_ext(G,10))

## Segundo Ejemplo: obtenido del artículo " M.E.J. Newman, Modularity and community structure 
## in networks, Proc. Natl. Acad. Sci. 103 (23) (2006) 8577–8582"
G1=nx.karate_club_graph()
mapping={}
N=G1.order()
for i in range(N):
    mapping[i]=i+1
G2=nx.relabel_nodes(G1,mapping)
## Represento gráficamente el grafo original
plt.figure(2)
nx.draw_shell(G2,with_labels=True)
print("Ejemplo2")
## Imprimo por pantalla la bipartición y la modularidad.
_,comunidades=opt_ext(G1,10)
colores=["green","lightgreen","red","lightcoral","brown","white"]
colors=[0 for i in range(N)]
for k in range(len(comunidades)):
    for i in comunidades[k]:
        colors[i]=colores[k]
plt.figure(figsize=(8, 6))
#biparticion,colores=normalized_min_cut(G2)
nx.draw_shell(G2,with_labels=True,node_color=colors,node_size=700,font_size=18)
plt.savefig('zachary_modularity_ext.jpg')
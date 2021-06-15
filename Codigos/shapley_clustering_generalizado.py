# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 23:15:10 2021

@author: 34625
"""
import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def shapley_community_detection(W,nodos,nc):
    try:
        # Tomo como variables los parámetros nodos, W_{ij}, nodos y k
        # nodos= lista de los nodos del grafo
        # k= número máximo de clusters en los que puede estar un nodo
        
        n=len(nodos)
        
        m = gp.Model("modelo")
        mini=min([min([W[i,j] for i in nodos if(W[i,j]>0)]) for j in nodos])
        
        # Introduzco las variables al modelo
        x={}
        z={}
        u={}
        for s in range(1,nc+1):
            for i in nodos:
                x[i,s]=m.addVar(vtype=GRB.BINARY, name="x"+str(i)+"_"+str(s))
        for s in range(1,nc+1):
            for i in nodos:
                for j in nodos:
                    if(i!=j):
                        z[i,j,s]=m.addVar(vtype=GRB.BINARY, name="z"+str(i)+"_"+str(j)+"_"+str(s))
                for r in range(s+1,nc+1):
                    u[i,s,r]=m.addVar(vtype=GRB.BINARY, name="u"+str(i)+"_"+str(s)+"_"+str(r))
        
        # Función objetivo
        m.setObjective(sum([sum([sum([W[i,j]*z[i,j,s] for j in nodos if(j!=i)]) for i in nodos]) for s in range(1,nc+1)]), GRB.MINIMIZE)
        
        # Restricciones
        for i in nodos:
            m.addConstr(sum([x[i,s] for s in range(1,nc+1)])>=1)
            m.addConstr(sum([x[i,s] for s in range(1,nc+1)])<=2)
        for s in range(1,nc):
            m.addConstr(sum([x[i,s] for i in nodos])>=sum([x[i,s+1] for i in nodos]))
        for i in nodos:
            for s in range(1,nc+1):
                m.addConstr(sum([z[i,j,s]*W[i,j] for j in nodos if(j!=i)])-sum([(x[i,s]-z[i,j,s])*W[i,j] for j in nodos if(j!=i)])>=0)
                for j in nodos:
                    if(i!=j):
                        m.addConstr(z[i,j,s] <= x[i,s])
                        m.addConstr(z[i,j,s] <= x[j,s])
                        m.addConstr(x[i,s]+x[j,s]-z[i,j,s]<=1)
                for r in range(s+1,nc+1):
                    m.addConstr(u[i,s,r] <= 1- x[i,s])
                    m.addConstr(u[i,s,r] <= x[i,r])
                    m.addConstr(x[i,r]-x[i,s]-u[i,s,r] <= 0)
        for s in range(1,nc+1):
            for r in range(s+1,nc+1):
                m.addConstr(n*sum([u[i,s,r] for i in nodos])>=sum([x[i,r] for i in nodos]))
        
        # Ejecuto el modelo
        m.optimize()
        
        # Imprimo las variables x_{is} que definen las comunidades
        for i in range(n*n):
            v=m.getVars()[i]
            if(v.x==1):
                print('%s %g' % (v.varName, v.x))
        
        comunidades=[]
        for s in range(1,nc+1):
            com=set()
            for i in nodos:
                if(x[i,s].x==1):
                    com.add(i)
            if(com!=set()):
                comunidades.append(com)
        print(comunidades)
        # Imprimo la funcíon objetivo
        print('Obj: %g' % m.objVal)
        color=["white" for i in range(n)]
        for j in comunidades[0]:
            color[nodos.index(j)]="red"
        for j in comunidades[1]:
            color[nodos.index(j)]="blue"
        for j in comunidades[2]:
            color[nodos.index(j)]="green"
        return color
    
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')

# Construyo el grafo a través de las aristas recopiladas en un txt
G6=nx.karate_club_graph()
mapping={}
N=G6.order()
for i in range(N):
    mapping[i]=i+1
grafo=nx.relabel_nodes(G6,mapping)
grafo=grafo.subgraph([1,2,3,4,5,6,7,9,10,15,16,33,32,24,25])
nodos=list(grafo.nodes())
n=len(nodos)
n_aristas=grafo.number_of_edges()
# grados de cada vértice
grados={}
for i in nodos:
    grados[i]=int(grafo.degree(i))
W1={}
for i in nodos:
    for j in nodos:
        if(i in grafo.neighbors(j)):
            W1[i,j]=1-grados[i]*grados[j]/(2*n_aristas)
        else:
            W1[i,j]=-grados[i]*grados[j]/(2*n_aristas)
    
# Construyo los parámetro P_ij y CN_ij
P={}
CN={}
for i in nodos:
    for j in nodos:
        if(i==j or grados[i]==0 or grados[j]==0):
            P[i,j]=0
            CN[i,j]=0
        else:
            P[i,j]=1/grados[i]+1/grados[j]
            vecinos_en_comun=set(grafo.neighbors(i)) & set(grafo.neighbors(j))
            CN[i,j]=(len(vecinos_en_comun)+1)*(1/grados[i]+1/grados[j])

# Construyo los parámetros W_ij
W={}
W_mod={}
for i in nodos:
    for j in nodos:
        if(i==j or grados[i]==0 or grados[j]==0):
            W[i,j]=0
        else:
            if(i not in grafo.neighbors(j)):
                W[i,j]=(CN[i,j]-P[i,j])/4
            else:
                if(grados[i]==1 or grados[j]==1):
                    W[i,j]=P[i,j]
                else:
                    W[i,j]=2*CN[i,j]+P[i,j]

# Ejecuto la función para k=3
color=shapley_community_detection(W1,nodos,6)
plt.figure(figsize=(8, 6))
nx.draw_shell(grafo,node_color=color,with_labels=True,node_size=700,font_size=18)
#plt.savefig('zachary_shapley_nc_6_k_2_com_6.jpg')
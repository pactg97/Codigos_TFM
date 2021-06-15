# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 17:05:59 2021

@author: 34625
"""

import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def modularity_overlapping(W,nodos,lam,nc):
    try:
        # Tomo como variables los parámetros nodos, W_{ij}, nodos y k
        # nodos= lista de los nodos del grafo
        # k= número máximo de clusters en los que puede estar un nodo
        
        n=len(nodos)
        
        m = gp.Model("modelo")
        #mini=min([min([W[i,j] for i in nodos if(W[i,j]>0)]) for j in nodos])
        W_total=sum([sum([W[i,j] for j in range(1,n+1)]) for i in range(1,n+1)])
        
        # Introduzco las variables al modelo
        x={}
        z={}
        u={}
        t={}
        h={}
        l1={}
        l2={}
        l3={}
        l4={}
        for s in range(1,nc+1):
            for i in range(1,n+1):
                x[i,s]=m.addVar(vtype=GRB.BINARY, name="x"+str(i)+"_"+str(s))
        for s in range(1,nc+1):
            for i in range(1,n+1):
                u[i,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="u"+str(i)+"_"+str(s))
                for r in range(s+1,nc+1):
                    h[i,s,r]=m.addVar(vtype=GRB.BINARY, name="h"+str(i)+"_"+str(s)+"_"+str(r))
        for s in range(1,nc+1):
            for i in range(1,n+1):
                for j in range(1,n+1):
                    z[i,j,s]=m.addVar(vtype=GRB.BINARY, name="z"+str(i)+"_"+str(j)+"_"+str(s))
#                    t[i,j,s]=m.addVar(vtype=GRB.BINARY, name="t"+str(i)+"_"+str(j)+"_"+str(s))
                    l1[i,j,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="l1"+"_"+str(i)+"_"+str(j)+"_"+str(s))
                    l2[i,j,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="l2"+"_"+str(i)+"_"+str(j)+"_"+str(s))
#                    l3[i,j,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="l3"+"_"+str(i)+"_"+str(j)+"_"+str(s))
#                    l4[i,j,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="l4"+"_"+str(i)+"_"+str(j)+"_"+str(s))
    
        # Función objetivo
        m.setObjective(sum([sum([sum([(l1[i,j,s]+l2[i,j,s])/2*(W[i,j]-sum([W[i,k] for k in range(1,n+1)])*sum([W[j,k] for k in range(1,n+1)])/W_total) for j in range(1,n+1)]) for i in range(1,n+1)]) for s in range(1,nc+1)]), GRB.MAXIMIZE)
        
        # Restricciones
            
        for s in range(1,nc):
            m.addConstr(sum([x[i,s] for i in range(1,n+1)])>=sum([x[i,s+1] for i in range(1,n+1)]))
        for i in range(1,n+1):
            for s in range(1,nc+1):
                m.addConstr(x[i,s]>=u[i,s]-lam)
                m.addConstr(x[i,s]<=1+u[i,s]-lam)
                for r in range(s+1,nc+1):
                    m.addConstr(sum([x[i,s] for s in range(1,nc+1)])>=1)
                    m.addConstr(sum([u[i,s] for s in range(1,nc+1)])==1)
                    m.addConstr(h[i,s,r] <= 1- x[i,s])
                    m.addConstr(h[i,s,r] <= x[i,r])
                    m.addConstr(x[i,r]-x[i,s]-h[i,s,r] <= 0)
                for j in range(1,n+1):
                    m.addConstr(z[i,j,s] <= x[i,s])
#                    m.addConstr(t[i,j,s] <= x[i,s])
                    m.addConstr(z[i,j,s] <= x[j,s])
#                    m.addConstr(t[i,j,s] <= 1-x[j,s])
                    m.addConstr(x[i,s]+x[j,s]-z[i,j,s]<=1)
#                    m.addConstr(x[i,s]-x[j,s]-t[i,j,s]<=0)
                    m.addConstr(l1[i,j,s] <= z[i,j,s])
                    m.addConstr(l2[i,j,s] <= z[i,j,s])
                    m.addConstr(l1[i,j,s] <= u[i,s])
                    m.addConstr(l2[i,j,s] <= u[j,s])
                    m.addConstr(l1[i,j,s] >= u[i,s]+z[i,j,s]-1)
                    m.addConstr(l2[i,j,s] >= u[j,s]+z[i,j,s]-1)
#                    m.addConstr(l3[i,j,s] <= x[i,s])
#                    m.addConstr(l3[i,j,s] <= 1-x[j,s])
#                    m.addConstr(l3[i,j,s] <= u[i,s])
#                    m.addConstr(l3[i,j,s] >= u[i,s]+x[i,s]-x[j,s]-1)
#                    m.addConstr(l4[i,j,s] <= x[i,s])
#                    m.addConstr(l4[i,j,s] <= 1-x[j,s])
#                    m.addConstr(l4[i,j,s] <= (1-u[j,s]))
#                    m.addConstr(l4[i,j,s] >= -u[j,s]+x[i,s]-x[j,s])
        for s in range(1,nc+1):
            for r in range(s+1,nc+1):
                m.addConstr(n*sum([h[i,s,r] for i in range(1,n+1)])>=sum([x[i,r] for i in range(1,n+1)]))
        # Ejecuto el modelo
        m.optimize()
        
        
        
        # Imprimo las variables x_{is} que definen las comunidades
        for i in range(n*nc):
            v=m.getVars()[i]
            if(v.x==1):
                print('%s %g' % (v.varName, v.x))
        comunidades=[]
        for s in range(1,nc+1):
            com=set()
            for i in range(1,n+1):
                print(u[i,s].varName," ",u[i,s].x)
                if(x[i,s].x==1):
                    com.add(i)
            if(com!=set()):
                comunidades.append(com)
        print([[nodos[i-1] for i in list(j)] for j in comunidades])
        color=["white" for i in range(n)]
        for j in comunidades[2]:
            color[j-1]="yellow"
        
        # Imprimo la funcíon objetivo
        print('Obj: %g' % m.objVal)
        return color
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')

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
W={}
for i in range(1,n+1):
    for j in range(1,n+1):
        if(nodos[i-1] in grafo.neighbors(nodos[j-1])):
            W[i,j]=1
        else:
            W[i,j]=0
#W1={}
#for i in nodos:
#    for j in nodos:
#        if(i in grafo.neighbors(j)):
#            W1[i,j]=1
#        else:
#            W1[i,j]=0
color=modularity_overlapping(W,nodos,0.25,4)
            
com=[{1,2,3,4,5,6,7,9,10,15,16,19,21,32,24,25,26}]
nc=len(com)
W_total=grafo.number_of_edges()
#plt.figure(figsize=(5, 3.5))
#plt.figure(figsize=(8, 6))
#nx.draw_shell(grafo,with_labels=True,node_color=color,node_size=700,font_size=18)
#plt.savefig('zachary_mod_ov_k_5_4.jpg')

#print(sum([sum([sum([W1[i,j] for j in com[s]]) for i in com[s]])-(sum([sum([W1[i,j]for j in nodos]) for i in com[s]]))*(sum([sum([W1[i,j]for j in nodos]) for i in com[s]]))/(2*W_total) for s in range(nc)]))
#
#print(4*sum([sum([sum([W1[i,j] for j in com[s]]) for i in com[s]])-(sum([sum([W1[i,j] for j in nodos]) for i in com[s]]))*(sum([sum([W1[i,j] for j in nodos]) for i in com[s]]))/(2*W_total) for s in range(nc)]))
#print(nodos)
#[[2, 3, 4, 9, 10], [33, 9], [32, 33, 15, 16, 24, 25], [7, 1, 5, 6]]
# 16.2898
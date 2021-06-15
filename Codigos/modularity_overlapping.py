# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 15:47:29 2021

@author: 34625
"""

import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def modularity_overlapping1(W,nodos,lam,nc):
    try:
        # Tomo como variables los parámetros nodos, W_{ij}, nodos y k
        # nodos= lista de los nodos del grafo
        # k= número máximo de clusters en los que puede estar un nodo
        
        n=len(nodos)
        
        m = gp.Model("modelo")
        mini=min([min([W[i,j] for i in range(1,n+1) if(W[i,j]>0)]) for j in range(1,n+1)])
        W_total=sum([sum([W[i,j] for j in range(1,n+1)]) for i in range(1,n+1)])
        
        # Introduzco las variables al modelo
        x={}
        z={}
        u={}
        l={}
        for s in range(1,n+1):
            for i in range(1,n+1):
                x[i,s]=m.addVar(vtype=GRB.BINARY, name="x"+str(i)+"_"+str(s))
        for s in range(1,nc+1):
            for i in range(1,n+1):
                u[i,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="u"+str(i)+"_"+str(s))
        for s in range(1,nc+1):
            for i in range(1,n+1):
                for j in range(1,n+1):
                    z[i,j,s]=m.addVar(vtype=GRB.BINARY, name="z"+str(i)+"_"+str(j)+"_"+str(s))
                    for k in range(1,n+1):
                        l[i,j,k,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="l"+str(i)+"_"+str(j)+"_"+str(k)+"_"+str(s))
    
        # Función objetivo
        m.setObjective(sum([sum([sum([W[i,j]*(u[i,s]+u[j,s])/2-sum([W[i,k]*(u[i,s]+u[k,s])/2 for k in range(1,n+1)])*sum([W[j,k]*(u[j,s]+u[k,s])/2 for k in range(1,n+1)])/(W_total) for j in range(1,n+1)]) for i in range(1,n+1)]) for s in range(1,nc+1)]), GRB.MAXIMIZE)
        
        # Restricciones
        for i in range(1,n+1):
            m.addConstr(sum([x[i,s] for s in range(1,nc+1)])>=1)
            m.addConstr(sum([u[i,s] for s in range(1,nc+1)])==1)
        for i in range(1,n+1):
            for s in range(1,nc+1):
                m.addConstr(x[i,s]>=u[i,s]-lam)
                m.addConstr(x[i,s]<=1+u[i,s]-lam)
                for j in range(1,n+1):
                    m.addConstr(z[i,j,s] <= x[i,s])
                    m.addConstr(z[i,j,s] <= x[j,s])
                    m.addConstr(x[i,s]+x[j,s]-z[i,j,s]<=1)
#                    for k in range(1,n+1):
#                        m.addConstr(l[i,j,k,s] <= z[i,j,s])
#                        m.addConstr(l[i,j,k,s] <= u[k,s]+1-x[k,s])
#                        m.addConstr(l[i,j,k,s] <= 1-u[k,s]+x[k,s])
#                        m.addConstr(l[i,j,k,s] >= u[k,s]+z[i,j,s]+x[k,s]-2)
#                        m.addConstr(l[i,j,k,s] >= z[i,j,s]-x[k,s]-u[k,s])
        
        # Ejecuto el modelo
        m.params.NonConvex=2
        m.optimize()
        
        # Imprimo las variables x_{is} que definen las comunidades
        for i in range(n*n):
            v=m.getVars()[i]
            if(v.x==1):
                print('%s %g' % (v.varName, v.x))
        
        # Imprimo la funcíon objetivo
        print('Obj: %g' % m.objVal)
    
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
        
def modularity_overlapping(W,nodos,lam,k,nc):
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
                    t[i,j,s]=m.addVar(vtype=GRB.BINARY, name="t"+str(i)+"_"+str(j)+"_"+str(s))
                    l1[i,j,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="l1"+str(i)+"_"+str(j)+"_"+str(s))
                    l2[i,j,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="l2"+str(i)+"_"+str(j)+"_"+str(s))
                    l3[i,j,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="l3"+str(i)+"_"+str(j)+"_"+str(s))
                    l4[i,j,s]=m.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="l4"+str(i)+"_"+str(j)+"_"+str(s))
    
        # Función objetivo
        m.setObjective(sum([sum([sum([W[i,j]*(l1[i,j,s]+l2[i,j,s])/2 for j in range(1,n+1)]) for i in range(1,n+1)])-(sum([sum([W[i,j]*(l1[i,j,s]+l2[i,j,s])/2 for j in range(1,n+1)]) for i in range(1,n+1)])+sum([sum([W[i,j]*(l3[i,j,s]+l4[i,j,s])/2 for j in range(1,n+1)]) for i in range(1,n+1)]))*(sum([sum([W[i,j]*(l1[i,j,s]+l2[i,j,s])/2 for j in range(1,n+1)]) for i in range(1,n+1)])+sum([sum([W[i,j]*(l3[i,j,s]+l4[i,j,s])/2 for j in range(1,n+1)]) for i in range(1,n+1)]))/W_total for s in range(1,nc+1)]), GRB.MAXIMIZE)
        
        # Restricciones
        for i in range(1,n+1):
            m.addConstr(sum([x[i,s] for s in range(1,nc+1)])>=1)
            #m.addConstr(sum([x[i,s] for s in range(1,nc+1)])<=k)
            m.addConstr(sum([u[i,s] for s in range(1,nc+1)])==1)
        for s in range(1,nc):
            m.addConstr(sum([x[i,s] for i in range(1,n+1)])>=sum([x[i,s+1] for i in range(1,n+1)]))
        for i in range(1,n+1):
            for s in range(1,nc+1):
                m.addConstr(x[i,s]>=u[i,s]-lam)
                m.addConstr(x[i,s]<=1+u[i,s]-lam)
                for r in range(s+1,nc+1):
                    m.addConstr(h[i,s,r] <= 1- x[i,s])
                    m.addConstr(h[i,s,r] <= x[i,r])
                    m.addConstr(x[i,r]-x[i,s]-h[i,s,r] <= 0)
                for j in range(1,n+1):
                    m.addConstr(z[i,j,s] <= x[i,s])
                    m.addConstr(t[i,j,s] <= x[i,s])
                    m.addConstr(z[i,j,s] <= x[j,s])
                    m.addConstr(t[i,j,s] <= 1-x[j,s])
                    m.addConstr(x[i,s]+x[j,s]-z[i,j,s]<=1)
                    m.addConstr(x[i,s]-x[j,s]-t[i,j,s]<=0)
                    m.addConstr(l1[i,j,s] <= z[i,j,s])
                    m.addConstr(l2[i,j,s] <= z[i,j,s])
                    m.addConstr(l1[i,j,s] <= u[i,s])
                    m.addConstr(l2[i,j,s] <= u[j,s])
                    m.addConstr(l1[i,j,s] >= u[i,s]+z[i,j,s]-1)
                    m.addConstr(l2[i,j,s] >= u[j,s]+z[i,j,s]-1)
                    m.addConstr(l3[i,j,s] <= x[i,s])
                    m.addConstr(l3[i,j,s] <= 1-x[j,s])
                    m.addConstr(l3[i,j,s] <= u[i,s])
                    m.addConstr(l3[i,j,s] >= u[i,s]+x[i,s]-x[j,s]-1)
                    m.addConstr(l4[i,j,s] <= x[i,s])
                    m.addConstr(l4[i,j,s] <= 1-x[j,s])
                    m.addConstr(l4[i,j,s] <= (1-u[j,s]))
                    m.addConstr(l4[i,j,s] >= -u[j,s]+x[i,s]-x[j,s])
        for s in range(1,nc+1):
            for r in range(s+1,nc+1):
                m.addConstr(n*sum([h[i,s,r] for i in range(1,n+1)])>=sum([x[i,r] for i in range(1,n+1)]))
        # Ejecuto el modelo
        m.optimize()
        
        
        # Imprimo las variables x_{is} que definen las comunidades
        for i in range(n*nc):
            v=m.getVars()[i]
            print('%s %g' % (v.varName, v.x))
#        for v in m.getVars():
#            print('%s %g' % (v.varName, v.x))
        
        # Imprimo la funcíon objetivo
        print('Obj: %g' % m.objVal)
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
grafo=grafo.subgraph([1,2,3,4,5,6,32,33,8,9])
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
W1={}
for i in nodos:
    for j in nodos:
        if(i in grafo.neighbors(j)):
            W1[i,j]=1
        else:
            W1[i,j]=0
modularity_overlapping(W,nodos,0.25,2,6)
            
com=[{32, 33, 9,1, 5, 6,8, 2, 3, 4}]
nc=len(com)
W_total=grafo.number_of_edges()
#print(sum([sum([sum([W1[i,j] for j in com[s]]) for i in com[s]])-(sum([sum([W1[i,j]for j in nodos]) for i in com[s]]))*(sum([sum([W1[i,j]for j in nodos]) for i in com[s]]))/(2*W_total) for s in range(nc)]))
#
#print(4*sum([sum([sum([W1[i,j]*0.25 for j in com[s]]) for i in com[s]])-(sum([sum([W1[i,j]*0.25 for j in nodos]) for i in com[s]]))*(sum([sum([W1[i,j]*0.25for j in nodos]) for i in com[s]]))/(2*W_total) for s in range(nc)]))
#print(nodos)
#print(4*sum([sum([sum([W[i,j]*0.25-sum([W[i,k]*0.25 for k in range(1,n+1)])*sum([W[j,k]*0.25 for k in range(1,n+1)])/(2*W_total) for j in range(1,n+1)]) for i in range(1,n+1)]) for s in range(1,nc+1)]))
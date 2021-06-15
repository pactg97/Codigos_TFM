# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 10:37:03 2020

@author: 34625
"""


import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def preproceso(G):
    n=G.order()
    mapping=dict()
    nodos=list(G.nodes())
    for i in range(n):
        mapping.update({nodos[i]:i})
    H= nx.relabel_nodes(G, mapping)
    return H

def calculo_B(A,n,m,k):
    B=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            B[i,j]=A[i,j]-k[i]*k[j]/(2*m)
    return B


def modularity_LP(G):
    try:
        # Create a new model
        m = gp.Model("mip1")
        n=G.order()
        n_aristas=G.number_of_edges()
        grados={}
        B={}
        nodos=list(G.nodes())
        aristas=list(G.edges())
        for (i,j) in G.degree():
            grados[i]=j
        for i in nodos:
            for j in nodos:
                if(i in G.neighbors(j)):
                    B[i,j]=1-grados[i]*grados[j]/(2*n_aristas)
                else:
                    B[i,j]=-grados[i]*grados[j]/(2*n_aristas)
        x={}
        
        for i in nodos:
            for j in nodos:
                x[i,j]=m.addVar(vtype=GRB.BINARY, name="x"+str(i)+"_"+str(j))
        m.setObjective(sum([sum([B[i,j]*x[i,j] for j in nodos if(i!=j)]) for i in nodos]), GRB.MAXIMIZE)
        
        for i in nodos:
            m.addConstr(x[i,i]==1)
            for j in nodos:
                m.addConstr(x[i,j]==x[j,i])
                for k in nodos:
                    if(i!=j and j!=k):
                        m.addConstr(x[i,j]+x[j,k]-x[i,k]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
                        m.addConstr(x[i,j]-x[j,k]+x[i,k]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
                        m.addConstr(-x[i,j]+x[j,k]+x[i,k]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
        # Optimize model
        m.optimize()
        print(len(m.getVars()))
        variables=m.getVars()
        k=0
        comunidades=[]
        for v in m.getVars():
            print('%s %g' % (v.varName, v.x))
        print('Obj: %g' % m.objVal)
        for i in nodos:
            conj={i}
            for j in nodos:
                if(variables[k].x==1):
                    conj.add(j)
                k=k+1
            comunidades.append(conj)
        comunidades_final=[]
        for i in comunidades:
            maximal=True 
            for j in [k for k in comunidades if k!=i]:
                if(i.issubset(j)):
                    maximal=False
            if(maximal):
                comunidades_final.append(i)
        
        com=[]
        for item in comunidades_final:
            if item not in com:
                com.append(item)
        return com,m.objVal
    
        
    
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
        
G1=nx.karate_club_graph()
mapping={}
N=G1.order()
for i in range(N):
    mapping[i]=i+1
G2=nx.relabel_nodes(G1,mapping)
grafo=G2.subgraph([1,2,3,4,5,6,7,9,10,15,16,33,32,24,25])
#comunidades,_=modularity_LP(G1)
#print(comunidades)
#colores=["green","lightgreen","red","lightcoral","brown","white"]
#colors=[0 for i in range(N)]
#for k in range(len(comunidades)):
#    for i in comunidades[k]:
#        colors[i]=colores[k]
#print(colors)
#
#plt.figure(figsize=(8, 6))
##biparticion,colores=normalized_min_cut(G2)
nx.draw_shell(grafo,with_labels=True)
##plt.savefig('zachary_modularity.jpg')

print(modularity_LP(grafo))
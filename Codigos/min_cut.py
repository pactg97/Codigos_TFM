# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 15:46:22 2021

@author: 34625
"""


import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import gurobipy as gp
from gurobipy import GRB

def calculo_B(A,n,m,k):
    B=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            B[i,j]=A[i,j]-k[i]*k[j]/(2*m)
    return B


def min_cut(G):
    try:
        A=nx.adjacency_matrix(G)
        # Create a new model
        m = gp.Model("mip1")
        n=G.order()
        n_aristas=G.number_of_edges()
        k=[G.degree(i) for i in range(n)]
        z=[]
        for i in range(n):
            z.append(m.addVar(vtype=GRB.BINARY, name="z"+str(i)))
        m.setObjective(sum([sum([A[i,j]*z[i]*(1-z[j]) for j in range(n)]) for i in range(n)]), GRB.MINIMIZE)
        m.addConstr(sum(z)>=1,"constraint")
        m.addConstr(sum(z)<=n-1,"constraint")
        # Optimize model
        m.optimize()
        print(len(m.getVars()))
        nx.draw_shell(G,with_labels=True)
        colores=[]
        for v in m.getVars():
            print('%s %g' % (v.varName, v.x))
            if(v.x==0):
                colores.append("green")
            else:
                colores.append("red")
        print('Obj: %g' % m.objVal)
        return colores
        
    
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
G6=nx.karate_club_graph()
mapping={}
N=G6.order()
for i in range(N):
    mapping[i]=i+1
G7=nx.relabel_nodes(G6,mapping)
colores=min_cut(G6)
plt.figure(2)
plt.figure(figsize=(8, 6))
nx.draw_shell(G7,with_labels=True,node_color=colores,node_size=700,font_size=18)
plt.savefig("zachary_min_cut")
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 10:29:37 2020

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

def min_clique_partition(G):
    try:
        A=nx.adjacency_matrix(G)
        # Create a new model
        m = gp.Model("mip1")
        n=G.order()
        x=[]
        for i in range(n):
            x0=[]
            for j in range(i):
                x0.append(m.addVar(vtype=GRB.BINARY, name="x"+str(j)+"_"+str(i)))
            x.append(x0)

        m.setObjective(sum([sum([(A[j,i]*(1-x[i][j])) for j in range(i)]) for i in range(n)]), GRB.MINIMIZE)
        
        for i in range(n):
            for j in range(i):
                m.addConstr(x[i][j]<=A[j,i],"c"+str(j)+"_"+str(i))
        for i in range(n):
            for j in range(i):
                for k in range(j):
                    m.addConstr(x[i][j]+x[j][k]-x[i][k]<=1,"c"+str(k)+"_"+str(j)+"_"+str(i))
        for i in range(n):
            for j in range(i):
                for k in range(j):
                    m.addConstr(x[i][j]-x[j][k]+x[i][k]<=1,"c"+str(k)+"_"+str(j)+"_"+str(i))
        for i in range(n):
            for j in range(i):
                for k in range(j):
                    m.addConstr(-x[i][j]+x[j][k]+x[i][k]<=1,"c"+str(k)+"_"+str(j)+"_"+str(i))
    
        # Optimize model
        m.optimize()
        print(len(m.getVars()))
        variables=m.getVars()
        nx.draw_shell(G,with_labels=True)
        k=0
        comunidades=[]
        for v in m.getVars():
            if(v.x==1):
                print('%s %g' % (v.varName, v.x))
        print('Obj: %g' % m.objVal)
        for i in range(n):
            conj={i}
            for j in range(i):
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
        return comunidades_final,m.objVal
    
        
    
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
        
G1=nx.Graph()
G1.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])
H=preproceso(G1)

G4=nx.Graph()
G4.add_nodes_from(range(34))
G4.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4),(1,5),(2,5),(3,5),(4,5),(6,5),(5,7),(5,8),(5,9),(5,10),(5,11),(5,12),
           (7,12),(8,10),(9,12),(10,11),(10,12),(10,13),(12,13),(5,13),(13,14),(11,14),(12,14),(5,14),(12,15),(5,15),
           (10,14),(11,12),(14,16),(13,25),(12,17),(14,26),(14,24),(14,22),(14,18),(15,25),(5,18),(5,19),(19,20),
            (20,21),(20,22),(19,21),(20,22),(21,23),(22,23),(17,18),(19,24),(16,25),(17,25),(17,26),(18,25),(18,26),
           (19,25),(19,26),(22,25),(23,25),(23,26),(27,25),(27,26),(28,25),(28,26),(23,29),(29,25),(29,26),(30,25),
           (30,26),(29,31),(25,31),(32,25),(32,26),(33,25),(33,26),(24,25)])

c1=min_clique_partition(G4)
c2=min_clique_partition(H)
print(c1)
print(c2)
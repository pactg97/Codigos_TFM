# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 10:33:07 2020

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

def bipartition(G,c):
    try:
        comunidades=[]
        # Create a new model
        n_arista_total=G.number_of_edges()
        k=[j for (i,j) in G.degree()]
        for comunidad in c:
            comunidad=list(comunidad)
            m = gp.Model("mip1")
            dt=sum([k[i] for i in comunidad])
            n=len(comunidad)
            aristas=[(u,v) for (u,v) in list(G.edges()) if (u in comunidad and v in comunidad)]
            n_arista=len(aristas)
            grados=[k[i] for i in comunidad]
            x0=[]
            x1=[]
            for i in range(n_arista):
                x0.append(m.addVar(vtype=GRB.BINARY, name="x0_"+str(i)))
                x1.append(m.addVar(vtype=GRB.BINARY, name="x1_"+str(i)))
            y=[]
            for i in range(n):
                y.append(m.addVar(vtype=GRB.BINARY, name="y"+str(i)))
            m.setObjective((sum(x1)+sum(x0))/n_arista_total-sum([grados[i]*y[i] for i in range(n)])*sum([grados[i]*y[i] for i in range(n)])/(4*n_arista_total*n_arista_total)-(dt-sum([grados[i]*y[i] for i in range(n)]))*(dt-sum([grados[i]*y[i] for i in range(n)]))/(4*n_arista_total*n_arista_total), GRB.MAXIMIZE)
            for i in range(n_arista):
                u,v=aristas[i]
                m.addConstr(x0[i]+x1[i]<=1,"cx_"+str(i))
                m.addConstr(2*x1[i]<=y[comunidad.index(u)]+y[comunidad.index(v)])
                m.addConstr(2-2*x0[i]>=y[comunidad.index(u)]+y[comunidad.index(v)])
            # Optimize model
            m.optimize()
            comunidad1=set()
            comunidad2=set()
            for v in m.getVars():
                print('%s %g' % (v.varName, v.x))
            for i in range(n):
                if(y[i].x==1):
                    comunidad1.add(comunidad[i])
                else:
                    comunidad2.add(comunidad[i])
            print('Obj: %g' % m.objVal)
            if(comunidad1 == set()):
                comunidades.append(comunidad2)
            else:
                if(comunidad2 == set()):
                    comunidades.append(comunidad1)
                else:
                    comunidades.append(comunidad1)
                    comunidades.append(comunidad2)          
        return comunidades
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
        
def division_heuristica(G):
    c=[set(range(G.order()))]
    comunidades=bipartition(G,c)
    while(c!=comunidades):
        c=comunidades
        print(c)
        comunidades=bipartition(G,c)
    return comunidades

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
c1=division_heuristica(G4)
c2=division_heuristica(H)
c3=division_heuristica(G1)
print(c1)
print(c2)
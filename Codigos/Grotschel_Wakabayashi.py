# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 17:47:22 2021

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

def fun_obj(c,clusters):
    sum=0
    for cluster in clusters:
        cluster=list(cluster)
        cluster=sorted(cluster)
        k=len(cluster)
        for i in range(k):
            for j in range(i+1,k):
                sum=sum+c[cluster[j]][cluster[i]]
    return sum

def Grotschel_Wakabayashi(M):
    try:
        n,p=M.shape
        # Create a new model
        m = gp.Model("mip1")
        x=[]
        c=[]
        for j in range(n):
            x0=[]
            c0=[]
            for i in range(j):
                x0.append(m.addVar(vtype=GRB.BINARY, name="x"+str(i)+"_"+str(j)))
                c0.append(len([k for k in range(p) if(M[i,k]!=M[j,k] and M[i,k]!="*" and M[j,k]!="*")])-len([k for k in range(p) if(M[i,k]==M[j,k] and M[j,k]!="*")]))
            x.append(x0)
            c.append(c0)
        m.setObjective(sum([sum([c[j][i]*x[j][i] for i in range(j)]) for j in range(n)]), GRB.MINIMIZE)
        
        for k in range(n):
            for j in range(k):
                for i in range(j):
                    if(c[j][i]<=0 or c[k][j]<=0):
                        m.addConstr(x[j][i]+x[k][j]-x[k][i]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
        for k in range(n):
            for j in range(k):
                for i in range(j):
                    if(c[j][i]<=0 or c[k][i]<=0):
                        m.addConstr(x[j][i]-x[k][j]+x[k][i]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
        for k in range(n):
            for j in range(k):
                for i in range(j):
                    if(c[k][j]<=0 or c[k][i]<=0):
                        m.addConstr(-x[j][i]+x[k][j]+x[k][i]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
    
        # Optimize model
        m.optimize()
        variables=m.getVars()
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
        
M=np.array([[0,5,0,0,0,1,2,3,3,4,1,"*",2,0,0],[0,4,0,0,3,0,4,3,3,0,0,"*",4,2,0],
   [0,4,0,0,3,3,4,3,3,3,0,"*",4,2,0],[1,2,1,2,2,0,0,0,2,0,1,1,2,1,2],
   [0,1,1,1,2,2,1,0,1,4,1,1,2,0,1],[1,3,1,0,0,1,1,1,1,2,0,1,2,0,2],
   [0,2,1,2,2,2,1,1,1,0,1,1,1,0,1],[0,1,0,0,0,3,3,3,3,3,0,"*",4,1,0],
   [0,5,0,0,0,1,2,3,3,1,1,"*",2,0,0],[0,3,1,0,2,3,1,0,1,3,1,1,1,0,1],
   [0,3,1,0,2,3,0,0,1,0,1,1,4,0,1],[1,2,1,2,2,0,0,0,2,0,1,1,4,1,2],
   [1,2,1,3,1,1,1,1,2,0,0,0,0,0,"*"],[0,0,0,0,1,0,0,0,0,0,1,"*",4,0,3],
   [0,2,1,2,2,2,1,1,1,0,1,1,2,0,1],[1,2,1,3,1,1,1,1,0,0,0,0,0,0,"*"],
   [0,1,1,1,0,2,1,1,1,0,1,1,1,0,1],[0,4,0,0,3,3,4,3,3,3,0,"*",4,2,0],
   [1,1,1,2,2,0,0,0,2,0,1,1,4,1,2],[1,3,1,0,0,1,0,0,1,3,0,1,2,0,2],
   [0,1,0,0,3,0,2,3,3,0,1,"*",2,0,0],[0,1,1,0,0,0,1,1,1,3,1,1,4,0,1],
   [1,3,1,0,2,0,1,1,1,1,1,1,3,0,1],[0,3,1,0,2,1,1,2,1,3,1,1,2,0,1],
   [0,1,1,0,1,0,1,1,1,0,1,1,4,0,1],[0,0,0,0,0,0,0,0,0,0,1,0,1,0,3],
   [1,2,1,3,1,1,1,1,2,1,0,0,0,0,1],[0,3,1,0,2,3,1,0,1,1,1,1,4,0,1],
   [0,2,1,2,2,2,1,1,1,2,1,1,3,0,1],[0,2,1,2,2,2,1,1,1,0,1,1,3,0,1],
   [0,2,1,2,2,2,1,1,1,3,1,1,4,0,1],[0,1,1,2,2,2,1,1,1,0,1,1,1,0,1],
   [1,2,1,3,1,1,1,1,2,1,0,0,0,0,"*"],[0,3,1,2,2,0,0,0,2,0,1,1,2,1,2],
   [0,2,1,2,2,2,1,1,1,0,1,1,4,0,1],[0,1,1,2,2,0,0,0,2,3,1,1,4,1,2]])
nombres=["Balaena","Balaenoptera","Balaenoptera Mus.","Berardius","Cephalorhynchus","Delphinapterus",
         "Delphinus","Eschrichtius","Eubalaena","Globicephala","Grampus","Hyperoodon",
         "Inia","Kogia","Lagenorhynchus","Lipotes","Lissodelphis","Megaptera","Mesoplodon",
         "Monodon","Neopbalaena","Neopbocaena","Orcaella","Orcinus","Phocaena","Physeter",
         "Platanista","Pseudorca","Sotalia","Sousa","Stenella","Steno","Stenodelphis",
         "Tasmacetus","Tursiops","Ziphius"]
comunidades,_=Grotschel_Wakabayashi(M)
print([{nombres[i] for i in comp} for comp in comunidades])



































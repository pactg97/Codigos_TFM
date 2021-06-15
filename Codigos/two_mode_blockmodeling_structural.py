# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 18:34:01 2021

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

def two_mode_blockmodeling_structural(S,B,P):
    try:
        N1,N2=S.shape
        K1,K2=B.shape
        # Create a new model
        m = gp.Model("mip1")
        v={}
        w={}
        eta={}
        z={}
        for k in range(K1):
            for i in range(N1):
                v[str(i)+","+str(k)]=m.addVar(vtype=GRB.BINARY, name="v"+str(i)+"_"+str(k))
        for k in range(K2):
            for i in range(N2):
                w[str(i)+","+str(k)]=m.addVar(vtype=GRB.BINARY, name="w"+str(i)+"_"+str(k))
        for i in range(N1):
            for j in range(N2):
                for k in range(K1):
                    for l in range(K2):
                        eta[str(i)+","+str(j)+","+str(k)+","+str(l)]=P[k,l]*(S[i,j]+B[k,l]-2*S[i,j]*B[k,l])
                        z[str(i)+","+str(j)+","+str(k)+","+str(l)]=m.addVar(vtype=GRB.CONTINUOUS,ub=1,lb=0, name="z"+str(i)+"_"+str(j)+"_"+str(k)+"_"+str(l))
        m.setObjective(sum([sum([sum([sum([eta[str(i)+","+str(j)+","+str(k)+","+str(l)]*z[str(i)+","+str(j)+","+str(k)+","+str(l)] for i in range(N1)]) for j in range(N2)]) for k in range(K1)])for l in range(K2)]), GRB.MINIMIZE)
        for i in range(N1):
            m.addConstr(sum([v[str(i)+","+str(k)] for k in range(K1)])==1,"1c_disjuntos_"+str(i))
        for k in range(K1):
            m.addConstr(sum([v[str(i)+","+str(k)] for i in range(N1)])>=1,"1c_novacios_"+str(k))
        for i in range(N2):
            m.addConstr(sum([w[str(i)+","+str(k)] for k in range(K2)])==1,"2c_disjuntos_"+str(i))
        for k in range(K2):
            m.addConstr(sum([w[str(i)+","+str(k)] for i in range(N2)])>=1,"2c_novacios_"+str(k))
        for i in range(N1):
            for j in range(N2):
                for k in range(K1):
                    for l in range(K2):
                        
                        m.addConstr(v[str(i)+","+str(k)]+w[str(j)+","+str(l)]-z[str(i)+","+str(j)+","+str(k)+","+str(l)]<=1,"c"+str(i)+","+str(j)+","+str(k)+","+str(l))
    
        # Optimize model
        m.optimize()
        particion1=[]
        particion2=[]
        for k in range(K1):
            particion1.append(set())
        for k in range(K2):
            particion2.append(set())
        variables=m.getVars()
        for v in variables:
            if(v.x==1):
                print('%s %g' % (v.varName, v.x))
        print('Obj: %g' % m.objVal)   
        l=0
        for k in range(K1):
            for i in range(N1):
                if(variables[l].x==1):
                    particion1[k].add(i)
                l=l+1
        for k in range(K2):
            for i in range(N2):
                if(variables[l].x==1):
                    particion2[k].add(i)
                l=l+1
        
        return particion1,particion2
            
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
        
southern_women=np.array([[1,1,1,1,1,1,0,1,1,0,0,0,0,0],
                         [1,1,1,0,1,1,1,1,0,0,0,0,0,0],
                         [0,1,1,1,1,1,1,1,1,0,0,0,0,0],
                         [1,0,1,1,1,1,1,1,0,0,0,0,0,0],
                         [0,0,1,1,1,0,1,0,0,0,0,0,0,0],
                         [0,0,1,0,1,1,0,1,0,0,0,0,0,0],
                         [0,0,0,0,1,1,1,1,0,0,0,0,0,0],
                         [0,0,0,0,0,1,0,1,1,0,0,0,0,0],
                         [0,0,0,0,1,0,1,1,1,0,0,0,0,0],
                         [0,0,0,0,0,0,1,1,1,0,0,1,0,0],
                         [0,0,0,0,0,0,0,1,1,1,0,1,0,0],
                         [0,0,0,0,0,0,0,1,1,1,0,1,1,1],
                         [0,0,0,0,0,0,1,1,1,1,0,1,1,1],
                         [0,0,0,0,0,1,1,0,1,1,1,1,1,1],
                         [0,0,0,0,0,0,1,1,0,1,1,1,0,0],
                         [0,0,0,0,0,0,0,1,1,0,0,0,0,0],
                         [0,0,0,0,0,0,0,0,1,0,1,0,0,0],
                         [0,0,0,0,0,0,0,0,1,0,1,0,0,0]])
supreme_court=np.array([[1,1,1,1,1,0,0,0,0],
                        [1,1,1,1,1,0,0,0,0],
                        [1,1,1,1,1,0,0,0,0],
                        [1,1,1,1,1,0,0,0,0],
                        [1,1,1,1,1,0,0,0,0],
                        [1,1,1,1,1,1,0,0,0],
                        [1,1,1,1,1,0,0,0,0],
                        [1,1,1,1,1,0,0,0,0],
                        [1,1,1,1,0,0,0,0,1],
                        [1,1,1,1,0,0,0,1,1],
                        [1,1,1,1,0,0,0,1,0],
                        [1,1,1,1,1,0,1,1,0],
                        [1,1,1,1,1,0,1,1,1],
                        [1,1,1,1,1,1,1,1,1],
                        [0,1,1,0,0,1,1,1,0],
                        [0,0,0,1,0,1,1,1,1],
                        [0,0,0,1,0,1,1,1,1],
                        [0,0,0,1,0,1,1,1,1],
                        [0,0,0,1,1,1,1,1,1],
                        [1,0,0,1,1,1,1,1,1],
                        [0,0,0,1,1,1,1,1,1],
                        [0,0,0,1,1,1,1,1,1],
                        [0,0,0,1,1,1,1,1,1],
                        [0,0,0,0,1,1,1,1,1],
                        [0,0,0,0,1,1,1,1,1],
                        [0,0,0,0,1,1,1,1,1]])
baker=open('Datos/baker.txt').read()
baker = [item.split() for item in baker.split('\n')[:-1]]
baker=np.array([[int(j) for j in i] for i in baker])
baker2=open('Datos/baker2.txt').read()
baker2 = [item.split() for item in baker2.split('\n')[:-1]]
baker2=np.array([[int(j) for j in i] for i in baker2])
B1=np.array([[1,1,0],[0,1,1]])
P1=np.array([[1,1,1],[1,1,1]])
B2=np.array([[1,0,0],[0,0,1]])
P2=P1
B3=np.array([[1,1,0],[0,1,0]])
P3=P1
B4=np.array([[1,1,0],[0,1,1],[0,1,0]])
P4=np.array([[1,1,100],[100,1,1],[100,1,100]])
B5=np.array([[0,1,1],[1,1,1],[1,1,0],[1,1,0],
             [1,1,0]])
P5=np.array([[1,1,1],[1,1,1],[1,1,1],[1,1,1],[1,1,1]])
B6=np.array([[0,1,1,1],[1,1,1,1],[1,1,1,0],[1,0,1,0],
             [1,1,0,0]])
P6=np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]])
B7=np.array([[0,1,1,1],[0,0,1,1],[1,1,1,1],[1,0,0,1],
             [1,1,1,0],[1,0,1,0],[1,1,0,0]])
P7=np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1],
             [1,1,1,1],[1,1,1,1],[1,1,1,1]])
B8=np.array([[0,1,0],[1,1,0]])
P8=np.array([[1,1,1],[1,1,1]])
B9=np.array([[0,0,0],[1,1,0],[0,1,0]])
P9=np.array([[1,1,1],[1,1,1],[1,1,1]])
B10=np.array([[0,1,0],[1,1,0],[0,1,1]])
P10=P9
indices_fila=[17,16,12,8,11,7,5,3,4,6,9,1,2,13,10,14,15,18]
indice_columna=[14,16,11,12,7,8,13,5,3,4,6,9,1,18,10,2,15,17]
ordenacion_columnas=[indice_columna.index(i) for i in indices_fila]

nueva_matriz=np.array([[baker2[i,j] for j in ordenacion_columnas] for i in range(18)])
print(nueva_matriz)
part1,part2=two_mode_blockmodeling_structural(nueva_matriz,B9,P9)
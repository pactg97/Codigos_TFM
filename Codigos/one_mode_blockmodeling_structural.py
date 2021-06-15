# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 10:09:47 2021

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

def one_mode_blockmodeling_structural(S,B,P):
    try:
        N=len(S)
        K=len(B)
        # Create a new model
        m = gp.Model("mip1")
        x={}
        delta={}
        y={}
        for i in range(N):
            for k in range(K):
                x[str(i)+","+str(k)]=m.addVar(vtype=GRB.BINARY, name="x"+str(i)+"_"+str(k))
        for i in range(N):
            for j in range(N):
                for k in range(K):
                    for l in range(K):
                        delta[str(i)+","+str(j)+","+str(k)+","+str(l)]=P[k,l]*(S[i,j]+B[k,l]-2*S[i,j]*B[k,l])
                        y[str(i)+","+str(j)+","+str(k)+","+str(l)]=m.addVar(vtype=GRB.CONTINUOUS,ub=1,lb=0, name="y"+str(i)+"_"+str(j)+"_"+str(k)+"_"+str(l))
        m.setObjective(sum([sum([sum([sum([delta[str(i)+","+str(j)+","+str(k)+","+str(l)]*y[str(i)+","+str(j)+","+str(k)+","+str(l)] for i in range(N)]) for j in range(N)]) for k in range(K)])for l in range(K)]), GRB.MINIMIZE)
        for i in range(N):
            m.addConstr(sum([x[str(i)+","+str(k)] for k in range(K)])==1,"c_disjuntos_"+str(i))
        for k in range(K):
            m.addConstr(sum([x[str(i)+","+str(k)] for i in range(N)])>=1,"c_novacios_"+str(k))
        for i in range(N):
            for j in range(N):
                for k in range(K):
                    for l in range(K):
                        m.addConstr(x[str(i)+","+str(k)]+x[str(j)+","+str(l)]-y[str(i)+","+str(j)+","+str(k)+","+str(l)]<=1,"c"+str(i)+","+str(j)+","+str(k)+","+str(l))
                        m.addConstr(y[str(i)+","+str(j)+","+str(k)+","+str(l)]<=x[str(i)+","+str(k)],"c"+str(i)+","+str(j)+","+str(k)+","+str(l))
                        m.addConstr(y[str(i)+","+str(j)+","+str(k)+","+str(l)]<=x[str(j)+","+str(l)],"c"+str(i)+","+str(j)+","+str(k)+","+str(l))
    
        # Optimize model
        m.optimize()
        particion=[]
        for k in range(K):
            particion.append(set())
        variables=m.getVars()
        s=100
        for i in range(N):
            for j in range(N):
                for k in range(K):
                    for l in range(K):
                        if(variables[s].x==1):
                            print('%s %g' % (variables[s].varName, variables[s].x),delta[str(i)+","+str(j)+","+str(k)+","+str(l)])
                        s=s+1

        print('Obj: %g' % m.objVal)   
        l=0
        for i in range(N):
            for k in range(K):
                if(variables[l].x==1):
                    particion[k].add(i)
                l=l+1
        
        return particion
            
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
        
transatlantic= np.array([[0,0,1,1,1,0,0,0,0,0,0,0,0],
  [1,0,1,0,0,0,0,0,0,0,1,0,0],
  [1,0,0,1,0,0,0,0,0,0,1,0,0],
  [1,1,1,0,0,0,0,0,0,0,0,0,0],
  [1,0,1,1,0,0,0,0,0,0,0,0,0],
  [0,0,0,0,1,0,0,0,0,0,0,1,1],
  [0,1,0,0,0,0,0,1,1,0,0,0,0],
  [0,1,0,0,0,0,1,0,1,0,0,0,0],
  [0,1,0,0,0,0,1,1,0,0,0,0,0],
  [1,0,0,0,0,1,0,0,0,0,0,0,0],
  [1,1,0,0,0,0,0,0,0,1,0,0,0],
  [1,0,0,0,0,1,0,0,0,1,0,0,0],
  [0,0,0,0,1,1,0,0,0,0,0,0,0]])

kansas_SAR= np.array([[0,1,1,1,0,1,1,0,1,0,1,0,0,0,1,1,1,0,0,1],
                      [1,0,0,1,1,1,1,0,0,0,1,1,1,1,0,1,1,0,0,1],
                      [1,1,0,1,0,1,1,0,1,0,0,0,0,1,0,1,1,0,0,0],
                      [1,1,0,0,1,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0],
                      [1,1,1,1,0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0],
                      [1,1,0,1,1,0,1,0,1,0,1,1,0,1,0,0,0,0,0,0],
                      [1,1,0,1,1,1,0,0,1,1,0,0,0,0,0,1,0,0,0,0],
                      [1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
                      [1,1,0,1,0,1,1,0,0,0,1,0,0,0,0,0,1,0,0,0],
                      [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [1,1,0,1,0,1,1,0,1,0,0,1,0,0,0,0,0,0,0,0],
                      [1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0],
                      [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                      [1,1,0,1,0,1,1,0,1,0,0,1,1,0,0,0,1,1,0,1],
                      [1,1,0,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,0,1],
                      [1,1,1,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0],
                      [1,1,1,1,0,1,1,0,1,1,0,0,0,1,0,0,0,0,0,0],
                      [1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0],
                      [1,1,1,0,0,1,1,0,1,0,0,0,0,0,0,0,1,0,0,0],
                      [1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0]])
## En el archivo .txt es necesario dejar una línea en blanco debajo de la matriz
## sino no lee la ultima fila de la matriz.
sharpstone=open('Datos/Sharpstone.txt').read()
sharpstone = [item.split() for item in sharpstone.split('\n')[:-1]]
sharpstone=np.array([[int(j) for j in i] for i in sharpstone])
polactor=open('Datos/PolActor.txt').read()
polactor = [item.split() for item in polactor.split('\n')[:-1]]
polactor=np.array([[int(j) for j in i] for i in polactor])
baker=open('Datos/baker.txt').read()
baker = [item.split() for item in baker.split('\n')[:-1]]
baker=np.array([[int(j) for j in i] for i in baker])
B1=np.array([[1,0,0,0],
             [0,0,0,0],
             [0,0,0,0],
             [0,1,0,1]])
P1=np.array([[1,1,1,1],
             [1,1,1,1],
             [1,1,1,1],
             [1,1,1,1]])
B2=np.array([[1,0,0,0],
             [1,0,0,0],
             [1,0,0,0],
             [0,1,0,1]])
B3=np.array([[1,0,0],
             [0,1,0],
             [0,0,0]])
P2=np.array([[1,1,1],
             [1,1,1],
             [1,1,1]])
B4=np.array([[1,1,1,0,0],
             [1,1,0,0,0],
             [1,1,0,0,0],
             [1,0,0,0,0],
             [1,1,1,1,0]])
P3=np.array([[1,1,1,1,1],
             [1,1,1,1,1],
             [1,1,1,1,1],
             [1,1,1,1,1],
             [1,1,1,1,1]])
P4=np.array([[1,1,1],[100,100,1],[1,100,1]])
kansas=open('Datos/kansas_sar.txt').read()
kansas = [item.split() for item in kansas.split('\n')[:-1]]
kansas=np.array([[int(j) for j in i] for i in kansas])
comunidades=one_mode_blockmodeling_structural(kansas_SAR,B4,P3)
nombres=["OSheriff","Highway P","Civil Def","Coroner","Attorney","Parks n Res","Game n Fish","Kansas DOT","Army Corps","Army Reserve","Crable Amb","Frank Co. Amb","Lee Rescue","Shawney", "Burl Police","Lynd Police","Red Cross", "Topica FD", "Carb FD","Topeka RBW"]
reordenacion=[]
for i in comunidades:
    reordenacion=reordenacion+list(i)
print([[kansas_SAR[i,j] for j in reordenacion] for i in reordenacion])
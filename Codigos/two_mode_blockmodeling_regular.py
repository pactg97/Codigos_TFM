# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 20:47:50 2021

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

def two_mode_blockmodeling_regular(S,B,P,lam,omega):
    try:
        N1,N2=S.shape
        K1,K2=B.shape
        # Create a new model
        m = gp.Model("mip1")
        v={}
        w={}
        alpha={}
        beta={}
        eta={}
        z={}
        for i in range(N1):
            for k in range(K1):
                v[str(i)+","+str(k)]=m.addVar(vtype=GRB.BINARY, name="v"+str(i)+"_"+str(k))
        for i in range(N2):
            for k in range(K2):
                w[str(i)+","+str(k)]=m.addVar(vtype=GRB.BINARY, name="w"+str(i)+"_"+str(k))
        for i in range(N1):
            for j in range(N2):
                for k in range(K1):
                    for l in range(K2):
                        eta[str(i)+","+str(j)+","+str(k)+","+str(l)]=P[k,l]*(S[i,j]+B[k,l]-2*S[i,j]*B[k,l])
                        z[str(i)+","+str(j)+","+str(k)+","+str(l)]=m.addVar(vtype=GRB.BINARY, name="z"+str(i)+"_"+str(j)+"_"+str(k)+"_"+str(l))
        for i in range(N1):
            for k in range(K1):
                for l in range(K2):
                    alpha[str(i)+","+str(k)+","+str(l)]=m.addVar(vtype=GRB.BINARY, name="alpha"+str(i)+"_"+str(k)+"_"+str(l))
        for i in range(N2):
            for k in range(K1):
                for l in range(K2):
                    beta[str(i)+","+str(k)+","+str(l)]=m.addVar(vtype=GRB.BINARY, name="beta"+str(i)+"_"+str(k)+"_"+str(l))
        m.setObjective(sum([sum([sum([sum([eta[str(i)+","+str(j)+","+str(k)+","+str(l)]*z[str(i)+","+str(j)+","+str(k)+","+str(l)] for i in range(N1) if(B[k,l]==0)]) for j in range(N2)]) for k in range(K1)])for l in range(K2)])+sum([sum([sum([alpha[str(i)+","+str(k)+","+str(l)]*lam[k,l] for i in range(N1) if(B[k,l]==1)]) for k in range(K1)])for l in range(K2)])+sum([sum([sum([beta[str(i)+","+str(k)+","+str(l)]*omega[k,l] for i in range(N2) if(B[k,l]==1)]) for k in range(K1)])for l in range(K2)]), GRB.MINIMIZE)
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
                            m.addConstr(v[str(i)+","+str(k)]+w[str(j)+","+str(l)]-z[str(i)+","+str(j)+","+str(k)+","+str(l)]<=1,"1c"+str(i)+","+str(j)+","+str(k)+","+str(l)) 
                            m.addConstr(v[str(i)+","+str(k)]+w[str(j)+","+str(l)]-2*z[str(i)+","+str(j)+","+str(k)+","+str(l)]>=0,"2c"+str(i)+","+str(j)+","+str(k)+","+str(l))
        for i in range(N1):
            for k in range(K1):
                for l in range(K2):
                    m.addConstr(sum([z[str(i)+","+str(j)+","+str(k)+","+str(l)]*S[i,j] for j in range(N2)])+alpha[str(i)+","+str(k)+","+str(l)]>=v[str(i)+","+str(k)],"3c"+str(i)+","+str(k)+","+str(l))
        for j in range(N2):
            for k in range(K1):
                for l in range(K2):
                    m.addConstr(sum([z[str(i)+","+str(j)+","+str(k)+","+str(l)]*S[i,j] for i in range(N1)])+beta[str(j)+","+str(k)+","+str(l)]>=w[str(j)+","+str(l)],"4c"+str(j)+","+str(k)+","+str(l))
        #m.addConstr(sum([sum([sum([sum([eta[str(i)+","+str(j)+","+str(k)+","+str(l)]*z[str(i)+","+str(j)+","+str(k)+","+str(l)] for i in range(N1) if(B[k,l]==0)]) for j in range(N2)]) for k in range(K1)])for l in range(K2)])+sum([sum([sum([alpha[str(i)+","+str(k)+","+str(l)]*lam[k,l] for i in range(N1) if(B[k,l]==1)]) for k in range(K1)])for l in range(K2)])+sum([sum([sum([beta[str(i)+","+str(k)+","+str(l)]*omega[k,l] for i in range(N2) if(B[k,l]==1)]) for k in range(K1)])for l in range(K2)])>=2)

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
        for i in range(N1):
            for k in range(K1):
                if(variables[l].x==1):
                    particion1[k].add(i)
                l=l+1
        for i in range(N2):
            for k in range(K2):
                if(variables[l].x==1):
                    particion2[k].add(i)
                l=l+1
        
        return particion1,particion2
            
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
        
baker=open('Datos/baker.txt').read()
baker = [item.split() for item in baker.split('\n')[:-1]]
baker=np.array([[int(j) for j in i] for i in baker])
baker2=open('Datos/baker2.txt').read()
baker2 = [item.split() for item in baker2.split('\n')[:-1]]
baker2=np.array([[int(j) for j in i] for i in baker2])
B1=np.array([[0,1,0],
             [1,1,0]])
P1=np.array([[1,100,1],
             [100,100,1]])
lam1=P1
omega1=P1
B2=np.array([[0,0,0],
             [1,1,0],
             [0,1,0]])
P2=np.array([[1,1,1],
             [100,100,1],
             [1,100,1]])
lam2=P2
omega2=P2
B3=np.array([[0,1,0],
             [1,1,0],
             [0,1,1]])
P3=np.array([[1,100,1],
             [100,100,1],
             [1,100,100]])
lam3=P3
omega3=P3
indices_fila=[17,16,12,8,11,7,5,3,4,6,9,1,2,13,10,14,15,18]
indice_columna=[14,16,11,12,7,8,13,5,3,4,6,9,1,18,10,2,15,17]
ordenacion_columnas=[indice_columna.index(i) for i in indices_fila]

nueva_matriz=np.array([[baker2[i,j] for j in ordenacion_columnas] for i in range(18)])

part1,part2=two_mode_blockmodeling_regular(nueva_matriz,B3,P3,P3,P3)
reordenacion1=[]
reordenacion2=[]
for k in range(3):
    reordenacion1=reordenacion1+list(part1[k])
    reordenacion2=reordenacion2+list(part2[k])
print(part1,part2)
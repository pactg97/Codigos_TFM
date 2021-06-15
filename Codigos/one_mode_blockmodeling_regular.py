# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 16:59:52 2021

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

def one_mode_blockmodeling_regular(S,B,P,lam,omega):
    try:
        N=len(S)
        K=len(B)
        # Create a new model
        m = gp.Model("mip1")
        x={}
        delta={}
        y={}
        alpha={}
        beta={}
        for k in range(K):
            for i in range(N):
                x[str(i)+","+str(k)]=m.addVar(vtype=GRB.BINARY, name="x"+str(i)+"_"+str(k))
        for i in range(N):
            for j in range(N):
                for k in range(K):
                    for l in range(K):
                        delta[str(i)+","+str(j)+","+str(k)+","+str(l)]=P[k,l]*(S[i,j]+B[k,l]-2*S[i,j]*B[k,l])
                        y[str(i)+","+str(j)+","+str(k)+","+str(l)]=m.addVar(vtype=GRB.BINARY, name="y"+str(i)+"_"+str(j)+"_"+str(k)+"_"+str(l))    
        for i in range(N):
            for k in range(K):
                for l in range(K):
                    if(B[k,l]==1):
                        alpha[str(i)+","+str(k)+","+str(l)]=m.addVar(vtype=GRB.BINARY, name="alpha"+str(i)+"_"+str(k)+"_"+str(l))
                        beta[str(i)+","+str(k)+","+str(l)]=m.addVar(vtype=GRB.BINARY, name="beta"+str(i)+"_"+str(k)+"_"+str(l))
        m.setObjective(sum([sum([sum([sum([delta[str(i)+","+str(j)+","+str(k)+","+str(l)]*y[str(i)+","+str(j)+","+str(k)+","+str(l)] for i in range(N)]) for j in range(N)]) for k in range(K) if(B[k,l]==0)])for l in range(K)])+sum([sum([sum([alpha[str(i)+","+str(k)+","+str(l)]*lam[k,l] for i in range(N)]) for k in range(K) if(B[k,l]==1)])for l in range(K)])+sum([sum([sum([beta[str(j)+","+str(k)+","+str(l)]*omega[k,l] for j in range(N)]) for k in range(K) if(B[k,l]==1)])for l in range(K)]), GRB.MINIMIZE)
        for i in range(N):
            m.addConstr(sum([x[str(i)+","+str(k)] for k in range(K)])==1,"c_disjuntos_"+str(i))
        for k in range(K):
            m.addConstr(sum([x[str(i)+","+str(k)] for i in range(N)])>=1,"c_novacios_"+str(k))   
        ## Tiene varias soluciones mínimas, para obtener la misma que el articulo
        ## fijo los tamaños de los cluster y fijo que la letra e (elemento 4) va al tercer cluster
        for i in range(N):
            for j in range(N):
                for k in range(K):
                    for l in range(K):
                        if(i!=j or k!=l):
                            m.addConstr(x[str(i)+","+str(k)]+x[str(j)+","+str(l)]-y[str(i)+","+str(j)+","+str(k)+","+str(l)]<=1,"1c"+str(i)+","+str(j)+","+str(k)+","+str(l))
                            m.addConstr(x[str(i)+","+str(k)]+x[str(j)+","+str(l)]-2*y[str(i)+","+str(j)+","+str(k)+","+str(l)]>=0,"2c"+str(i)+","+str(j)+","+str(k)+","+str(l))
        for i in range(N):
            for k in range(K):
                for l in range(K):
                    if(B[k,l]==1):
                        m.addConstr(sum([y[str(i)+","+str(j)+","+str(k)+","+str(l)]*S[i,j] for j in range(N)])+alpha[str(i)+","+str(k)+","+str(l)]>=x[str(i)+","+str(k)],"3c"+str(i)+","+str(k)+","+str(l))
        for j in range(N):
            for k in range(K):
                for l in range(K):
                    if(B[k,l]==1):
                        m.addConstr(sum([y[str(i)+","+str(j)+","+str(k)+","+str(l)]*S[i,j] for i in range(N)])+beta[str(j)+","+str(k)+","+str(l)]>=x[str(j)+","+str(l)],"4c"+str(j)+","+str(k)+","+str(l))
        # Optimize model
        m.optimize()
        particion=[]
        for k in range(K):
            particion.append(set())
        variables=m.getVars()
        for v in variables:
            if(v.x==1):
                print('%s %g' % (v.varName, v.x))
        print('Obj: %g' % m.objVal)   
        l=0
        for k in range(K):
            for i in range(N): 
                if(variables[l].x==1):
                    particion[k].add(i)
                l=l+1
        return particion  
            
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
        
        
everett=open('everett.txt').read()
everett = [item.split() for item in everett.split('\n')[:-1]]
everett=np.array([[int(j) for j in i] for i in everett])
baker=open('baker.txt').read()
baker = [item.split() for item in baker.split('\n')[:-1]]
baker=np.array([[int(j) for j in i] for i in baker])
baker2=open('baker2.txt').read()
baker2 = [item.split() for item in baker2.split('\n')[:-1]]
baker2=np.array([[int(j) for j in i] for i in baker2])
B1=np.array([[1,1,0],
             [1,0,1],
             [0,1,1]])
P1=np.array([[1,1,1],
             [1,1,1],
             [1,1,1]])
lam1=P1
omega1=P1
B2=np.array([[0,0,0],[1,1,0],[0,1,0]])
P2=np.array([[1,1,1],[100,100,1],[1,100,1]])
lam2=P2
omega2=lam2
B3=np.array([[0,1,0],[1,1,0],[0,1,1]])
P3=np.array([[1,100,1],[100,100,1],[1,100,100]])
lam3=P3
omega3=P3
print(everett)
indices_fila=[17,16,12,8,11,7,5,3,4,6,9,1,2,13,10,14,15,18]
indice_columna=[14,16,11,12,7,8,13,5,3,4,6,9,1,18,10,2,15,17]
ordenacion_columnas=[indice_columna.index(i) for i in indices_fila]

nueva_matriz=np.array([[baker2[i,j] for j in ordenacion_columnas] for i in range(18)])
comunidades=one_mode_blockmodeling_regular(nueva_matriz,B3,P3,P3,P3)
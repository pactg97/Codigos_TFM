# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 11:26:11 2021

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

def MTZ2_graph_connected_M(M,G):
    try:
        E=list(G.edges())
        n,p=M.shape
        # Create a new model
        m = gp.Model("mip1")
        x={}
        y={}
        c={}
        f={}
        l={}
        A=[(i,j) for (i,j) in E]+[(j,i) for (i,j) in E]
        for k in range(n):
            y[str(k)]=m.addVar(vtype=GRB.BINARY, name="z"+str(k))
        for j in range(n):
            for i in range(j):
                x[str(i)+","+str(j)]=m.addVar(vtype=GRB.BINARY, name="x"+str(i)+"_"+str(j))
                c[str(i)+","+str(j)]=len([k for k in range(p) if(M[i,k]!=M[j,k])])-len([k for k in range(p) if(M[i,k]==M[j,k])])
        for i,j in A:
            f[str(i)+","+str(j)]=m.addVar(vtype=GRB.BINARY ,name="f"+str(i)+"_"+str(j))
        for i in range(n):
            l[str(i)]=m.addVar(vtype=GRB.BINARY, name="l_"+str(i))
        m.setObjective(sum([sum([c[str(i)+","+str(j)]*x[str(i)+","+str(j)] for i in range(j)]) for j in range(n)]), GRB.MINIMIZE)
        for k in range(n):
            for j in range(k):
                for i in range(j):
                    if(c[str(i)+","+str(j)]<=0 or c[str(j)+","+str(k)]<=0):
                        m.addConstr(x[str(i)+","+str(j)]+x[str(j)+","+str(k)]-x[str(i)+","+str(k)]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
        for k in range(n):
            for j in range(k):
                for i in range(j):
                    if(c[str(i)+","+str(j)]<=0 or c[str(i)+","+str(k)]<=0):
                        m.addConstr(x[str(i)+","+str(j)]-x[str(j)+","+str(k)]+x[str(i)+","+str(k)]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
        for k in range(n):
            for j in range(k):
                for i in range(j):
                    if(c[str(j)+","+str(k)]<=0 or c[str(i)+","+str(k)]<=0):
                        m.addConstr(-x[str(i)+","+str(j)]+x[str(j)+","+str(k)]+x[str(i)+","+str(k)]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
        for (i,j) in A:
            m.addConstr(l[str(i)]+1<=l[str(j)]+n*(1-f[str(i)+","+str(j)]))
            if(i<j):
                m.addConstr(f[str(i)+","+str(j)]+f[str(j)+","+str(i)]<=x[str(i)+","+str(j)])
            else:
                m.addConstr(f[str(i)+","+str(j)]+f[str(j)+","+str(i)]<=x[str(j)+","+str(i)])
        for k in range(n):
            m.addConstr(y[str(k)]+sum([x[k,j] for j in range(k+1,n)])>=1)
            m.addConstr(sum([f[str(i)+","+str(k)] for i in range(n) if((i,k) in A)])==1-y[str(k)])
        for j in range(n):
            for k in range(j):
                m.addConstr(y[str(k)]+x[str(k)+","+str(j)])
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

def MTZ2_graph_connected_cost(C,G):
    try:
        E=list(G.edges())
        n=len(C)
        # Create a new model
        m = gp.Model("mip1")
        x={}
        y={}
        c={}
        f={}
        l={}
        A=E+[(j,i) for (i,j) in E]
        for k in range(n):
            y[str(k)]=m.addVar(vtype=GRB.BINARY, name="y"+str(k))
        for j in range(n):
            for i in range(j):
                x[str(i)+","+str(j)]=m.addVar(vtype=GRB.BINARY, name="x"+str(i)+"_"+str(j))
                c[str(i)+","+str(j)]=C[i,j]
        for i,j in A:
            f[str(i)+","+str(j)]=m.addVar(vtype=GRB.BINARY,name="f"+str(i)+"_"+str(j))
        for i in range(n):
            l[str(i)]=m.addVar(vtype=GRB.CONTINUOUS,lb=0, name="l_"+str(i))
        m.setObjective(sum([sum([c[str(i)+","+str(j)]*x[str(i)+","+str(j)] for i in range(j)]) for j in range(n)]), GRB.MINIMIZE)
        for k in range(n):
            for j in range(k):
                for i in range(j):
                    if(c[str(i)+","+str(j)]<=0 or c[str(j)+","+str(k)]<=0):
                        m.addConstr(x[str(i)+","+str(j)]+x[str(j)+","+str(k)]-x[str(i)+","+str(k)]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
        for k in range(n):
            for j in range(k):
                for i in range(j):
                    if(c[str(i)+","+str(j)]<=0 or c[str(i)+","+str(k)]<=0):
                        m.addConstr(x[str(i)+","+str(j)]-x[str(j)+","+str(k)]+x[str(i)+","+str(k)]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
        for k in range(n):
            for j in range(k):
                for i in range(j):
                    if(c[str(j)+","+str(k)]<=0 or c[str(i)+","+str(k)]<=0):
                        m.addConstr(-x[str(i)+","+str(j)]+x[str(j)+","+str(k)]+x[str(i)+","+str(k)]<=1,"c"+str(i)+"_"+str(j)+"_"+str(k))
        for (i,j) in A:
            m.addConstr(l[str(i)]+1<=l[str(j)]+n*(1-f[str(i)+","+str(j)]))
            if(i<j):
                m.addConstr(f[str(i)+","+str(j)]+f[str(j)+","+str(i)]<=x[str(i)+","+str(j)])
            else:
                m.addConstr(f[str(i)+","+str(j)]+f[str(j)+","+str(i)]<=x[str(j)+","+str(i)])
        for k in range(n):
            m.addConstr(y[str(k)]+sum([x[str(k)+","+str(j)] for j in range(k+1,n)])>=1)
            m.addConstr(sum([f[str(i)+","+str(k)] for i in range(n) if((i,k) in A)])==1-y[str(k)])
        for j in range(n):
            for k in range(j):
                m.addConstr(y[str(k)]+x[str(k)+","+str(j)]<=1)
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
        
M1=np.array([[0,5,0,0,0,1,2,3,3,4,1,"*",2,0,0],[0,4,0,0,3,0,4,3,3,0,0,"*",4,2,0],
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

G1=nx.complete_graph(36)
matriz_C=open('Datos/matriz_C.txt').read()
matriz_C = [item.split() for item in matriz_C.split('\n')[:-1]]
matriz_C=np.array([[int(float(j)) for j in i] for i in matriz_C])
G2=nx.complete_graph(20)
aristas=open('Datos/grafo.txt').read()
aristas= [item.split() for item in aristas.split('\n')[:-1]]
aristas=[tuple([int(i[0])-1,int(i[1])-1]) for i in aristas]
G3=nx.Graph()
G3.add_nodes_from(range(20))
G3.add_edges_from(aristas)
MTZ2_graph_connected_cost(matriz_C,G3)
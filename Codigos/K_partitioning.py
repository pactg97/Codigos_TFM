# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 10:26:47 2021

@author: 34625
"""
import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import gurobipy as gp
from gurobipy import GRB

matriz_C=open('Datos/distancias_baviera.txt').read()
matriz_C = [item.split() for item in matriz_C.split('\n')[:-1]]
matriz_C=np.array([[int(float(j)) for j in i] for i in matriz_C])

matriz=[[matriz_C[i,j] for i in range(20)] for j in range(20)]

def K_partitioning_cost(C,K):
    try:
        n=len(C)
        # Create a new model
        m = gp.Model("mip1")
        x={}
        z={}
        c={}
        for k in range(n):
            for i in range(k+1):
                z[i,k]=m.addVar(vtype=GRB.BINARY, name="z"+str(i)+"_"+str(k))
        for j in range(n):
            for i in range(j):
                x[i,j]=m.addVar(vtype=GRB.BINARY, name="x"+str(i)+"_"+str(j))
                c[i,j]=C[i][j]
        m.setObjective(sum([sum([c[i,j]*x[i,j] for i in range(j)]) for j in range(n)]), GRB.MINIMIZE)
        
        for i in range(n):
            for j in range(i+1,n):
                for k in range(i,n):
                    m.addConstr(x[i,j]<=z[i,k]+sum([z[j,l] for l in range(j,n) if(l!=k)]),"c"+str(i)+"_"+str(j)+"_"+str(k))
        for k in range(n):
            for j in range(k+1):
                for i in range(j):
                    m.addConstr(x[i,j]>=z[i,k]+z[j,k]-1,"1c"+str(i)+"_"+str(j)+"_"+str(k))
        for k in range(n):
            for i in range(k+1):
                m.addConstr(z[i,k]<=z[k,k])
        for i in range(n):
            m.addConstr(sum([z[i,k] for k in range(i,n)])==1)
        m.addConstr(sum([z[k,k] for k in range(n)])==K)
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
            conj=set()
            for j in range(i+1):
                if(variables[k].x==1):
                    conj.add(j+1)
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

print(K_partitioning_cost(matriz,4))
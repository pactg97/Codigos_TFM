# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 11:53:19 2021

@author: 34625
"""
import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import math
from sklearn.metrics import pairwise_distances
from sklearn.datasets import make_blobs
import matplotlib.pyplot as plt

n=500
d=2
jain,y=make_blobs(n, d, centers=2,random_state=10)

W=np.zeros([n,n])
for i in range(len(jain)):
    for j in range(len(jain)):
        W[i,j]=math.exp(-np.linalg.norm(jain[i,:]-jain[j,:])/5)

def transform(x,A):
    n=len(x)
    x_final=np.zeros(n)
    for i in range(n):
        x_final[i]=x[i]*(np.dot(A,x)[i])/(np.dot(np.transpose(x),np.dot(A,x)))
    return x_final



def dominant_set(A, x=None, epsilon=1.0e-10):
    """Compute the dominant set of the similarity matrix A with the
    replicator dynamics optimization approach. Convergence is reached
    when x changes less than epsilon.
    See: 'Dominant Sets and Pairwise Clustering', by Massimiliano
    Pavan and Marcello Pelillo, PAMI 2007.
    """
    if x is None:
        x = np.ones(A.shape[0])/float(A.shape[0])
        
    distance = epsilon*2
    while distance > 0:
        x_old = x.copy()
        # x = x * np.dot(A, x) # this works only for dense A
        x = x * A.dot(x) # this works both for dense and sparse A
        x = x / x.sum()
        distance = np.linalg.norm(x - x_old)

    return x

conj=dominant_set(W)
plt.figure(figsize=(8, 5))
plt.scatter(jain[:,0],jain[:,1])
datos=[]
for i in range(n):
    if(conj[i]>0):
        datos.append(jain[i,:])

datos=np.array(datos)
plt.scatter(datos[:,0],datos[:,1])
#plt.savefig('dominant_set2.jpg')





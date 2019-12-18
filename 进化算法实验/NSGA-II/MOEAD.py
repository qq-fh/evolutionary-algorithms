# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 16:09:58 2019

@author: Feng hui
"""

import numpy as np
import random
import GA
import NSGA_II 
import time


'''
功能：获得随机的权重向量
ObjVar：目标变量数量
N：目标问题数量
'''
def getRandomVector(ObjVar, N):
    weightsVec = []
    for i in range(N):
        weightsVec.append(np.random.rand(ObjVar))
        weightsVec[i] /= np.sum(weightsVec[i])
    return np.array(weightsVec)

#获得权重向量的邻居
def getNeighbor(weightsVec, neiNum):
    n = len(weightsVec)
    B = np.zeros((n,neiNum))
    dis = np.zeros((n,1))
    for i in range(n):
        for j in range(n):
            dis[j] = np.sum( (weightsVec[i]-weightsVec[j])**2 )
        sortIndex = np.argsort(dis, axis=0)   #升序排列
        B[i,:] = sortIndex[0:neiNum].T  #取前neiNum个近邻
    return B

#初始化解，并且初始化每个目标函数的最优解
def initOptiSol(boundList, func, solNums, encodeLength, encodeAccuracy):
    n = len(boundList)   #变量数量
    Sols = GA.getIntialPopulation(encodeLength, solNums)
    decodeSols = GA.decodedChromosome(encodeLength, Sols, boundList, encodeAccuracy)
    Objs = func(decodeSols)
    Z = np.max(Objs, axis=0)
    return Z, Sols
'''
获得初始化解(实数)
'''
def getActualInitSol(N , boundList):
    var = len(boundList)   #变量的数量
    Sols = np.zeros((N , var))
    for i in range(N):
        for j in range(var):
            Sols[i, j] = boundList[j][0]+np.random.rand()*(boundList[j][1]-boundList[j][0])
    return Sols
            
    
def weightsFunc(x, Z, weights):
    vec = np.abs((x-Z)*weights)
    if vec.ndim==1:
        return np.sum(vec)
    else:
        return np.sum(vec, axis=1)

#判断x是否被y支配
def isDominated(x, y):
    if (y<x).any() and (y<=x).all():
        return True
    else :
        return False
#输入两个向量  返回支配者
def returnDominator(x1, x2):
    if isDominated(x1, x2): #x1被x2支配
        return x2,1
    else:
        return x1,0
#直接交叉
def crossoverDirect(x):
    m,n = x.shape
    updatepopulation = x.copy()
    crossoverPoint = random.sample(range(1, n), 1)
    crossoverPoint = crossoverPoint[0]   #随机选择交叉点
    # one-single-point crossover
    updatepopulation[0, 0:crossoverPoint] = x[0, 0:crossoverPoint]
    updatepopulation[0, crossoverPoint: ] = x[1, crossoverPoint:]
    updatepopulation[1, 0:crossoverPoint] = x[1, 0:crossoverPoint]
    updatepopulation[1, crossoverPoint: ] = x[0, crossoverPoint:]
    return updatepopulation
'''
实数交叉
'''
def actualCrossover(p1, p2, boundList):
    x1 = p1.copy()
    x2 = p2.copy() 
    yita = 1   #交叉参数
    
    if (x1==x2).all():   #两个解相同没必要交叉
        return x1,x2
    
    for i in range(len(x1)):
        u = np.random.rand()

        u = np.random.rand()
        if( u<0.5 ):
            r = (2*u)**(1/(yita+1))
        else:
            r = (1/(2*(1-u)))**(1/(yita+1))
            
        x1[i] = 0.5*( (1+r)*x1[i] + (1-r)*x2[i] )
        x2[i] = 0.5*( (1-r)*x1[i] + (1+r)*x2[i] )
           
        x1[i] = min(x1[i], boundList[i][1])    #修复交叉后的解
        x1[i] = max(x1[i], boundList[i][0])
        x2[i] = min(x2[i], boundList[i][1])
        x2[i] = max(x2[i], boundList[i][0])
    return x1,x2

def actualCrossover2(population, boundList, Pc):
    m, n = population.shape
    numbers = np.uint8(m * Pc)
    # 确保进行交叉的染色体个数是偶数个
    if numbers % 2 != 0:
        numbers += 1
    # 交叉后得到的新种群
    updatepopulation = np.zeros((m, n), dtype=np.uint8)
    # 产生随机索引
    index = random.sample(range(m), numbers)
    # 不进行交叉的染色体进行复制
    for i in range(m):
        if not index.__contains__(i):
            updatepopulation[i, :] = population[i, :]
    # crossover
    while len(index) > 0:
        a = index.pop()
        b = index.pop()
        updatepopulation[a],updatepopulation[b] = actualCrossover(population[a], population[b], boundList)
    return updatepopulation
'''
多项式变异 实数
'''
def polyMutation(p, boundList):
    x = p.copy()
    yita = 1   #变异参数
    for i in range(len(x)):
        r = np.random.rand()
        if r>3/len(p):
#        if r>0.1:
            continue
        u = np.random.rand()
        if u<0.5:
            delta = (2*u)**(1/(yita+1))-1
        else:
            delta = 1 - ( (2*(1-u))**(1/(yita+1)) )
        x[i] += delta
        x[i] = min(x[i], boundList[i][1])  #修复变异解
        x[i] = max(x[i], boundList[i][0])
    return x
def polyMutation2(P, boundList):
    m,n = P.shape
    updatePopu = np.zeros((m,n))
    for i in range(m):
        updatePopu[i] = polyMutation(P[i], boundList)
    return updatePopu

def ZDT1(x):
    m,n = x.shape  #决策变量个数
    f1 = x[:,0].reshape(m,1)     
    g = 9*np.sum(x[:,1:], axis=1)/(n-1)
    g = g.reshape(m,1)+1
    f2 = g*(1-(f1/g)**0.5)
    return np.concatenate((f1,f2), axis=1)   


def ZDT2(x):
    m,n = x.shape  #决策变量个数
    f1 = x[:,0].reshape(m,1)     
    g = 9*np.sum(x[:,1:], axis=1)/(n-1)
    g = g.reshape(m,1)+1
    f2 = g*(1-(f1/g)**2)
    return np.concatenate((f1,f2), axis=1)   #必须是两维的数据才能聚合
    
def ZDT3(x):
    m,n = x.shape  #决策变量个数
    f1 = x[:,0].reshape(m,1)     
    g = 9*np.sum(x[:,1:], axis=1)/(n-1)
    g = g.reshape(m,1)+1
    f2 = g*(1-(f1/g)**0.5 - (f1/g)*np.sin(10*np.pi*f1))
    return np.concatenate((f1,f2), axis=1)   #必须是两维的数据才能聚合

def generate_2D_mean_vector(N):
    vec = np.zeros((N,2))
    step = 1/N
    for i in range(N):
        vec[i,0] = step*i
        vec[i,1] = 1-vec[i,0]
    return vec







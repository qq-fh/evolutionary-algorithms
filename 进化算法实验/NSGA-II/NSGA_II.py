# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 09:02:37 2019

@author: Administrator
"""

import numpy as np
import matplotlib.pyplot as plt 
import random
import time

#目标函数1
def func1(x):
    return x**2

#目标函数2
def func2(x):
    return (x-2)**2

#输入x为1维实数
#输出f为2维实数
def SCH(x):
    f1 = x**2
    f2 = (x-2)**2
    return np.concatenate((f1,f2), axis=1)
#输入x为3维实数
#输出f为2维实数
def FON(x):
    n = len(x)
    f = np.zeros((n,2))
    for i in range(n):
        f[i,0] = 1-np.exp(-np.sum( (x[i,:]-1/np.sqrt(3))**2) )
        f[i,1] = 1-np.exp(-np.sum( (x[i,:]+1/np.sqrt(3))**2) )
    return f
#输入x为2维
#输出f为2维
def KUR(x):
    n,var = x.shape
    f1 = np.zeros((n,1))
    f2 = np.zeros((n,1))
    for i in range(n):
        for j in range(var-1):
            f1[i] += -10*np.exp( -0.2*np.sqrt(x[i,j]**2 + x[i,j+1]**2) )
        f2[i] = np.sum(  np.abs(x[i,:])**0.8 + 5*np.sin(x[i,:]**3) )
    return np.concatenate( (f1,f2), axis=1)  

#快速非支配排序算法
#vals：目标函数值 （num维）
def fast_non_dominated_sort(vals):
    m,num = vals.shape  #获得目标函数值维数  样本数量
    S  = [[] for i in range(m)]  #被个体p支配的个体集合
    n = np.zeros((m, 1))            #支配p的个体数
    rank = np.zeros((m,1))   #每个个体的分层
    front = [[]]      #非支配前沿
    for p in range(m):  #对每一个个体
        S[p] = []
        n[p] = 0
        for q in range(m):
            pval = vals[p,:]  #取出 p q的值
            qval = vals[q,:]
            if((pval<qval).any() and (pval<=qval).all()):  #p中的分量全部小于等于q中的分量  而且至少有一个小于 
                if q not in S[p]:
                    S[p].append(q)   #q被p支配 
            elif ((pval>qval).any() and (pval>=qval).all()):
                n[p]+=1
        if n[p]==0:
            rank[p] = 0
            front[0].append(p)   #第一级前沿
    i = 0
    while (front[i] != []):   #当前最前沿不为空集
        Q = []
        for p in front[i]:   
            for q in S[p]:
                n[q] -= 1
                if n[q]==0:
                    rank[q] = i+1
                    if q not in Q:
                        Q.append(q)
        i += 1 
        front.append(Q)
    del front[-1]    #删除最后一个空集
    return front,rank


#计算某一个前沿的拥挤度  为实现经经营策略奠定基础
#I：为某一前沿
def crowding_distance_assignment(I, front):
    I = I[front]
    m,n = I.shape
    iDis = np.zeros((m,1))    #每个个体的拥挤度
    iDis[0]  = 44444444444444   #边界的拥挤度为无穷大
    iDis[-1] = 444444444444444
    for i in range(n):
        I = I[np.argsort(I[:, i])]   #升序排列 
        for j in range(1,m-1):
            iDis[j] += I[j+1, i]-I[j-1, i]
    return iDis
    
    
    

    
    
    
    
    
    
                



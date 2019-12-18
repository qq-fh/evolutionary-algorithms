# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:43:09 2019

@author: Administrator
"""

import numpy as np
from scipy.optimize import fsolve
import random


def getChromosomeLen(boudList, delta):
    lengths = []  #染色体的长度
    for i in boudList:
        lower = i[0]   #下界  
        upper = i[1]
        length = fsolve( lambda x:((upper-lower)/delta)-2**x+1, 50 )
        lengths.append( int(np.floor(length)) ) #取整
    return lengths


def decodeChromosomes(lengths, chromosones, boudList, delta = 0.0001 ):
    N = chromosones.shape[0]   #染色体数量
    varNum = len(lengths)   #变量数量
    decodeVals = np.zeros( (N , varNum) )  #解码之后的解
    for i,chromo in enumerate(chromosones):       #i为下标
        start = 0
        for index,L in enumerate(lengths):
            power = L-1  #指数权
            demical = 0
            for k in range(start, start+int(L)):
                demical += chromo[k]* (2**power)
                power-=1
            lower = boudList[index][0]  
            upper = boudList[index][1]
            decodeVal = lower + (upper-lower)*demical/(2**L-1)
            decodeVals[i,index] = decodeVal
            start += L
    return decodeVals    #返回解码后的染色体


'''
轮赌算法选择新的种群
'''
def SelNewPopu(chromosomes, fitNessSum):
    N,m = chromosomes.shape
    randoms = np.random.rand(N)   #随机生成N的0-1的小数
    newPopu = np.zeros( (N,m) )
    for i in range(N):
        logical = randoms[i]<=fitNessSum
        index = np.where(logical == 1)
        newPopu[i] = chromosomes[index[0][0], :]
    return newPopu
    
'''
单点交叉函数
'''
def crossover(chromosomes, pc=0.8):
    N,geNum = chromosomes.shape #染色体数量  基因数量
    crossNum = np.int(N*pc)   #交叉个体的数量
    crossNum = (crossNum+1)//2 * 2  #保证偶数个交叉
    newCrossPopu = np.zeros( (N, geNum) )  #交叉过后的种群
    index = random.sample(range(N), crossNum)  #生成交叉的下标
    for i in range(N):
        if i not in index:
            newCrossPopu[i,:] = chromosomes[i,:]   #没被选择的染色体直接复制
    while len(index)>0:
        a = index.pop()   #返回最顶上的值  并且删除此值
        b = index.pop()
        crossPoint = random.sample(range(geNum), 1)[0]   #随机生成交叉点
        newCrossPopu[a, 0:crossPoint] = chromosomes[a, 0:crossPoint]
        newCrossPopu[b, 0:crossPoint] = chromosomes[b, 0:crossPoint]
        newCrossPopu[a, crossPoint:] = chromosomes[b, crossPoint:]     #单点交叉
        newCrossPopu[b, crossPoint:] = chromosomes[a, crossPoint:]
    return newCrossPopu

def mutation(chromosomes, pm=0.01):
    m,n = chromosomes.shape
    mutationPopu = chromosomes.copy()   #复制原始染色体
    muGeneNum = np.int(pm*m*n)  #变异基因的数量
    index = random.sample(range(m*n), muGeneNum)  #随机生成变异基因下标
    for gene in index:
        chromIndex = gene//n   #染色体的位置
        geneIndex = gene%n     #基因的位置
        mutationPopu[chromIndex, geneIndex] = 1- mutationPopu[chromIndex, geneIndex]  #基因翻转变异
    return mutationPopu  #返回经过变异的染色体



def fitFunc(x):
    return 21.5 + x[0]*np.sin(4*np.pi*x[0]) + x[1]*np.sin(20*np.pi*x[1])
def fitFunc1(x):
    n = len(x)
    S = 0
    for i in range(n):
        S += x[i]*np.sin(np.sqrt(abs(x[i])))
    return S

def getFitness(func, decodeVals):   
    N = decodeVals.shape[0]     
    fitnessVals = np.zeros( (N ,1) )  #适应度
    for i in range(N):
        fitnessVals[i] = func(decodeVals[i])
    fitnessValsStd = fitnessVals/sum(fitnessVals)
    fitSum = np.zeros( (N,1) )
    fitSum [0] = fitnessValsStd[0]
    for i in range(1,N):
        fitSum[i] += fitSum[i-1]+fitnessValsStd[i] 
    return fitnessVals, fitSum








    
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
#        if(u>0.8):
#            continue
        
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
'''
多项式变异 实数
'''
def polyMutation(p, boundList):
    x = p.copy()
    yita = 1   #变异参数
    for i in range(len(x)):
        r = np.random.rand()
        if r>1/len(p):
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

def DTLZ1(x):
    m,n = x.shape
    g = 100*(n-2) + 100*np.sum( (x[:,2:]-0.5)**2-np.cos(20*np.pi*(x[:,2:]-0.5)) , axis=1)
    f1 = (1+g)*x[:,0]*x[:,1]
    f2 = (1+g)*x[:,0]*(1-x[:,1])
    f3 = (1+g)*(1-x[:,0])
    return np.concatenate((f1.reshape(m,1),f2.reshape(m,1),f3.reshape(m,1)), axis=1)

def DTLZ2(x):
    m,n = x.shape
    g = np.sum(x[:, 2:]**2, axis=1)
    f1 = (1+g)*np.cos(np.pi*x[:,0]/2)*np.cos(np.pi*x[:,1]/2)
    f2 = (1+g)*np.cos(np.pi*x[:,0]/2)*np.sin(np.pi*x[:,1]/2)
    f3 = (1+g)*np.sin(np.pi*x[:,0]/2)
    return np.concatenate((f1.reshape(m,1),f2.reshape(m,1),f3.reshape(m,1)), axis=1)


def generate_2D_mean_vector(N):
    vec = np.zeros((N,2))
    step = 1/N
    for i in range(N):
        vec[i,0] = step*i
        vec[i,1] = 1-vec[i,0]
    return vec

def MOEAD_Algorithm(Func, iterNum, nEP, N, T, isPrint):

    ObjVar = 2             #目标变量的个数
    codeAccuracy = 0.0001  #编码精度
    EP = []                #外部种群
    var = 30                #输入变量维数
    bound = [0,1]
    boundList = []
    for i in range(var):
        boundList.append(bound)
    encodeLengths = GA.getEncodedLength(codeAccuracy, boundList)
    
    nowTim = 0       
    sTim = time.time()    #记录程序开始运行的系统时间
    #---------------------------------------初始化-----------------------------------
    weightsVec = generate_2D_mean_vector(N)       #获取均匀分布的权重向量
    weightsNeighbor = getNeighbor(weightsVec, T)  #计算向量之间的近似度 并且返回近邻向量的下标
    Z, Sols = initOptiSol(boundList, Func, N, encodeLengths, codeAccuracy ) 
    #---------------------------------------更新--------------------------------------
    for iters in range(iterNum):  
        for i in range(N):   
            #----------------------随机在近邻中抽取两个个体--------------------
            rdIdx = random.sample(range(T), 2)      #随机采取两个近邻
            weightsNeighbor = weightsNeighbor.astype(np.uint16)   #修改数据类型为int
            randomNei = weightsNeighbor[i,rdIdx]     #array数据形式 可以使用非整数
            x = Sols[randomNei, :]                   #array数据类型可以像matlab一样操作
            
            #----------------------------交叉变异--------------------------------
            #单点交叉
            crossovers = crossoverDirect(x)    
            #位点变异
            mutations = GA.mutation(crossovers, Pm = 0.01)   
            #染色体解码
            decodePopu = GA.decodedChromosome(encodeLengths, mutations, boundList)   
            
            #-----------------------------更新Z--------------------------------
            #计算新解对应的目标函数值向量
            vals = Func(decodePopu)                       
            #返回支配者,对应的下标
            newSolVal, dominatorIdx = returnDominator(vals[0], vals[1])  
            Idx = newSolVal>Z
            Z[Idx] = newSolVal[Idx]   
            
            #----------------------------更新临近解------------------------------
            #当前个体的近邻解
            neighbors = Sols[weightsNeighbor[i], :]  
            #近邻染色体解码
            decodeNeighbors = GA.decodedChromosome(encodeLengths, neighbors, boundList)  
            #计算目标函数值
            ObjVals = Func(decodeNeighbors)   
            #计算新解对应的加权函数值
            newObj = weightsFunc(newSolVal, Z, weightsVec[i])
            #新解对应的染色体
            newSolChromo = mutations[dominatorIdx, :]                  
            #所有临近解随影的加权函数值
            weightsObjs  = np.zeros((len(neighbors), 1)) 
#            for k in range(len(neighbors)):
#                weightsObjs[k, 0] = weightsFunc(ObjVals[k, :], Z, weightsVec[weightsNeighbor[k]])
#            weightsObjs = weightsFunc(ObjVals, Z, weightsVec[i])
            
            for k in range(T):
                #判断产生的新个体是否比近邻要优秀
#                if weightsObjs[k]<=newObj:     
                if  weightsFunc(ObjVals[k, :], Z, weightsVec[weightsNeighbor[i,k]]) <= weightsFunc(newSolVal, Z, weightsVec[weightsNeighbor[i,k]]):
                    Sols[weightsNeighbor[i,k], :] = newSolChromo
    
            exist = False
            delete = []  #需要删除的解
            for k in range(len(EP)):
                if (newSolVal==EP[k]).all(): #解已经存在
                    exist = True
                    break
                if isDominated(EP[k], newSolVal):  #被新个体支配
                    delete.append(k)
                elif isDominated(newSolVal, EP[k]):  #新个体被支配
                    exist=True
            delete.reverse()
            if len(delete)>0:
                for k in delete:  #删除被支配的个体
                    del EP[k]
            if exist==False:
                EP.append(newSolVal)
            while len(EP)>nEP:
                selected = np.random.randint(0, len(EP))  #生成一个随机整数 删除EP中的多余解
                del EP[selected]
               
        #是否打印运行状态
        if(isPrint):   
            if (time.time()-nowTim)>=3:
                nowTim = time.time()
                schedule = 1-(iterNum-iters)/iterNum
                print('MOEAD--运行时间：%.3f s----进度：%.2f%%----EP：%d'%(nowTim-sTim, schedule*100, len(EP)))
    return EP, len(EP), time.time()-sTim








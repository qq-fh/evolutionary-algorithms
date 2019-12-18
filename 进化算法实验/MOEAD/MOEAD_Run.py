# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 16:11:13 2019

@author: feng hui
"""
import numpy as np
import matplotlib.pyplot as plt
import MOEAD
import time
import random
import GA

Func = MOEAD.ZDT3
#Func = NSGA_II.KUR
iterNum = 200        #迭代次数  
nEP = 200              #外部存档的数量
N = 100                #生成子问题的个数
T = 20                 #每个子问题的紧邻个数
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
sTim = time.time()
#---------------------------------------初始化-----------------------------------
#weightsVec = np.loadtxt('test.txt')
weightsVec = MOEAD.generate_2D_mean_vector(N)       #获取均匀分布的权重向量
weightsNeighbor = MOEAD.getNeighbor(weightsVec, T)  #计算向量之间的近似度 并且返回近邻向量的下标
Z, Sols = MOEAD.initOptiSol(boundList, Func, N, encodeLengths, codeAccuracy ) 
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
        crossovers = MOEAD.crossoverDirect(x)    
        #位点变异
        mutations = GA.mutation(crossovers, Pm = 0.01)   
        #染色体解码
        decodePopu = GA.decodedChromosome(encodeLengths, mutations, boundList)   
        
        #-----------------------------更新Z--------------------------------
        #计算新解对应的目标函数值向量
        vals = Func(decodePopu)                       
        #返回支配者,对应的下标
        newSolVal, dominatorIdx = MOEAD.returnDominator(vals[0], vals[1])  
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
        newObj = np.sum( weightsVec[weightsNeighbor[i], :]* np.abs(newSolVal-Z ) , axis=1)
        #新解对应的染色体
        newSolChromo = mutations[dominatorIdx, :]      
        
        #所有临近解对应的加权函数值
        weightsObjs = np.sum( weightsVec[weightsNeighbor[i] , :]*np.abs(ObjVals-Z), axis=1 )
      
        CompareIndex = newObj>weightsObjs
        Sols[weightsNeighbor[i, CompareIndex], :] = newSolChromo

        exist = False
        delete = []  #需要删除的解
        for k in range(len(EP)):
            if (newSolVal==EP[k]).all(): #解已经存在
                exist = True
                break
            if MOEAD.isDominated(EP[k], newSolVal):  #被新个体支配
                delete.append(k)
            elif MOEAD.isDominated(newSolVal, EP[k]):  #新个体被支配
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
           
           
    if (time.time()-nowTim)>=3:
        nowTim = time.time()
        schedule = 1-(iterNum-iters)/iterNum
        print('MOEAD--运行时间：%.3f s----进度：%.2f%%----EP：%d'%(nowTim-sTim, schedule*100, len(EP)))
    

#-------------------------绘图-------------------------
x1 = [x[0] for x in EP]
x2 = [x[1] for x in EP]
plt.scatter(x1, x2, color='r', label='PF' )
plt.xlabel('func1')
plt.ylabel('func2')
plt.title('MOEAD_ZDT3')
plt.legend()
plt.show()









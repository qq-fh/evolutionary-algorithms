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
from mpl_toolkits.mplot3d import Axes3D

#选择是否保存外部群体
SAVE_EP = 0  

Func = MOEAD.DTLZ2
iterNum = 200          #迭代次数  
nEP = 200              #外部存档的数量
N = 300                #生成子问题的个数
T = 20                 #每个子问题的紧邻个数
ObjVar = 3             #目标变量的个数
EP = []                #外部种群
var = 10               #输入变量维数
bound = [0,1]
boundList = []
if Func==MOEAD.DTLZ2:  #测试函数的边界不一样
    for i in range(2):
        boundList.append([0,1])
    for i in range(var-2):
        boundList.append([-1,1])
else:
    for i in range(var):
        boundList.append(bound)

nowTim = 0
sTim = time.time()
#---------------------------------------初始化-----------------------------------
weightsVec = MOEAD.getRandomVector(ObjVar, N)
#weightsVec = np.loadtxt('test.txt')
weightsNeighbor = MOEAD.getNeighbor(weightsVec, T)  #计算向量之间的近似度 并且返回近邻向量的下标
weightsNeighbor = weightsNeighbor.astype(np.uint16) #修改数据类型为int
Sols = MOEAD.getActualInitSol(N, boundList)         #初始化实数解
InitFVals = Func(Sols)                              #初始化目标函数值 
Z = np.min(InitFVals , axis=0)                      #初始化Z    理想最小点
#---------------------------------------更新--------------------------------------
for iters in range(iterNum):  
    for i in range(N):   
        #----------------------随机在近邻中抽取两个个体--------------------
        rdIdx = random.sample(range(T), 2)          #随机采取两个近邻
        randomNeighbor = weightsNeighbor[i,rdIdx]   #近邻下标
        x = Sols[randomNeighbor, :]                 #取出近邻解 
        
        #----------------------------交叉变异--------------------------------
        #实数交叉模拟二进制交叉  SBX
        crossovers = MOEAD.actualCrossover(x[0], x[1], boundList)   
        #多项式变异
        mutation1 = MOEAD.polyMutation(crossovers[0], boundList)  
        mutation2 = MOEAD.polyMutation(crossovers[1], boundList)
        mutation1 = mutation1.reshape(1, len(mutation1))
        mutation2 = mutation2.reshape(1, len(mutation2))

        #-----------------------------更新Z--------------------------------
        #计算新解对应的目标函数值向量     
        val1 = Func(mutation1)                       
        val2 = Func(mutation2)
        #返回支配者,对应的下标
        newSolVal, dominatorIdx = MOEAD.returnDominator(val1, val2)  
        newSolVal = newSolVal.reshape(newSolVal.shape[1],)
        Idx = newSolVal<Z
        Z[Idx] = newSolVal[Idx]   
        
        #----------------------------更新临近解------------------------------
        #当前个体的近邻解
        neighbors = Sols[weightsNeighbor[i], :]  
        #计算目标函数值
        ObjVals = Func(neighbors)   
        #计算新解对应的加权函数值
        newWeightObj = np.max( weightsVec[weightsNeighbor[i], :]* np.abs(newSolVal-Z ) , axis=1)
        #新解对应的染色体
        if dominatorIdx==0:
            newSol = mutation1
        else:
            newSol= mutation2      
        #所有临近解对应的加权函数值
        weightsObjs = np.max( weightsVec[weightsNeighbor[i] , :]*np.abs(ObjVals-Z), axis=1 )
      
        CompareIndex = newWeightObj < weightsObjs
        Sols[weightsNeighbor[i, CompareIndex], :] = newSol

        #----------------------------保存外部种群------------------------------
        if SAVE_EP:         #是否保存EP
            exist = False
            delete = []     #需要删除的解
            for k in range(len(EP)):
                if (newSolVal==EP[k]).all(): #解已经存在
                    exist = True
                    break
                if MOEAD.isDominated(EP[k], newSolVal):   #被新个体支配
                    delete.append(k)
                elif MOEAD.isDominated(newSolVal, EP[k]): #新个体被支配
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
           
    if (time.time()-nowTim)>=1:
        nowTim = time.time()
        schedule = 1-(iterNum-iters)/iterNum
        print('MOEAD--运行时间：%.3f s----进度：%.2f%%----EP：%d'%(nowTim-sTim, schedule*100, len(EP)))
        
    

#-------------------------绘图-------------------------
x1 = np.array([x[0] for x in EP]).reshape(len(EP), 1)
x2 = np.array([x[1] for x in EP]).reshape(len(EP), 1)
x3 = np.array([x[2] for x in EP]).reshape(len(EP), 1)

fig = plt.figure()
f1 = fig.add_subplot(121, projection='3d')

f1.scatter(x1, x2, x3, color='b',marker='o', label='EP' )
f1.set_xlabel('f1')
f1.set_ylabel('f2')
f1.set_zlabel('f3')

plt.title('MOEAD')
plt.legend()

Fvs = Func(Sols)
x1 = np.array([x[0] for x in Fvs]).reshape(len(Fvs), 1)
x2 = np.array([x[1] for x in Fvs]).reshape(len(Fvs), 1)
x3 = np.array([x[2] for x in Fvs]).reshape(len(Fvs), 1)

f2 = fig.add_subplot(122, projection='3d')
f2.scatter(x1, x2, x3, color='r', marker='o', label='FV' )
f2.set_xlabel('f1')
f2.set_ylabel('f2')
f2.set_zlabel('f3')
plt.title('MOEAD')
plt.legend()

plt.show()

#-------------------------保存数据------------------------
VALS = np.concatenate((x1,x2,x3), axis=1)
if   Func==MOEAD.ZDT1:
    np.savetxt('ZDT1.txt', VALS, fmt='%10.6f')
elif Func==MOEAD.ZDT2:
    np.savetxt('ZDT2.txt', VALS, fmt='%10.6f')
elif Func==MOEAD.ZDT3:
    np.savetxt('ZDT3.txt', VALS, fmt='%10.6f')
elif Func==MOEAD.DTLZ1:
    np.savetxt('DTLZ1.txt', VALS, fmt='%10.6f')
elif Func==MOEAD.DTLZ2:
    np.savetxt('DTLZ2.txt', VALS, fmt='%10.6f')



# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:24:11 2019

@author: Administrator
"""

import numpy as np
import matplotlib.pyplot as plt
import GA
import time

stim = time.time()


N = 10 #初始种群数量
n = 30 #变量维数
K = 3  #进行仿真的函数编号
bouds = [[-500,500], [-5.12,5.12], [-32,32], [-600,600], [-100,100]]  #
func = [GA.fitness_func1, GA.fitness_func2, GA.fitness_func3, GA.fitness_func4, GA.fitness_func9]
boudList = [] #自变量的边界
for i in range(n):
    boudList.append(bouds[K-1])
    
lengthsEncode = GA.getEncodedLength(boundarylist=boudList)  #获得编码的长度

optiFits = [] #每次迭代获得的最优适应度
optiSols = [] #最优的解

population = GA.getIntialPopulation(lengthsEncode, N)   #随机初始化种群

for iterition in range(20000):     #进行迭代进化
    chromDecode =  GA.decodedChromosome(lengthsEncode, population, boudList)  #解码染色体
    fitness, fitnessSum = GA.getFitnessValue(func[K-1], chromDecode)  #计算当前种群的适应度
    newSelPopu = GA.selectNewPopulation(population, fitnessSum)     #轮赌选择新的群体 1/3
    crossoverPopu = GA.crossover(newSelPopu, Pc = 0.5)              #单点交叉         2/3
    
#    chromDecode =  GA.decodedChromosome(lengthsEncode, crossoverPopu, boudList)  #解码染色体
#    fitnessMid, fitnessSumMid = GA.getFitnessValue(func[K-1], chromDecode)  #计算当前种群的适应度
#    mutationPopu = GA.AdaptiveMutation(crossoverPopu, fitnessMid, [0.001, 0.01])             #变异             3/3
    
    mutationPopu = GA.mutation(crossoverPopu, Pm = 0.01)
    
    chromDecode2 = GA.decodedChromosome(lengthsEncode, mutationPopu, boudList)  #解码
    fitness2, fitnessSum2 = GA.getFitnessValue(func[K-1], chromDecode2) #计算交叉变异选择之后的群体适应度
    
    fitnessMerge = np.concatenate((fitness, fitness2))  #混合操作之前跟之后的群体适应度  选出较优的个体
    chromMerge = np.concatenate( (population, mutationPopu) )
    selIndex = fitnessMerge.argsort(axis=0)   #适应度从小到大排列的下标
    selIndex = [np.int(x) for x in selIndex]
    population = chromMerge[selIndex[:N], :]   #新的群体诞生
    
    minFit = np.min(fitnessMerge)
    optiFits.append(minFit)
    print("iter :%d---optiFitness：%f---meanFitness：%f"%(iterition, minFit,np.mean(fitness)))

etim = time.time()
print("程序运行时间：%d s"%(etim-stim))
np.savetxt("optiFitnesss_func4.txt", optiFits,fmt="%10.4f")

plt.figure(1)
plt.plot(optiFits, c='b', label='Fitness with iteration')
plt.legend()
plt.show()



    
    
    
    
    
    
    
    
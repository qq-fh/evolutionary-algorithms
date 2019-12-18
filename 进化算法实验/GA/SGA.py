# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 10:21:09 2019

@author: Administrator
"""

import numpy as np
from scipy.optimize import fsolve
import GA_Func as ga
import matplotlib.pyplot as plt
import time


start= time.time()


'''
初始给出的条件
'''
delta = 0.0001    #求解的精度
#boudList = [ [-3.0, 12.1], [4.1, 5.8] ]   #解的定义域  边界
boudList = []   
N = 10      #初始群体数量
n = 10       #变量维数
for i in range(n):
    boudList.append([-32, 32])   #生成变量边界

lengths = ga.getChromosomeLen(boudList, delta)  #染色体的长度

chromosomes = np.zeros((N ,int(sum(lengths))))   #染色体初始化
for i in range(N):  #初始化种群
    chromosomes[i, :] = np.random.randint(0, 2, int(sum(lengths)))
    

optiSolutions = []   #优化得到的解
optiVals = []        #优化得到的适应度   可以看优化曲线
for it in range(5000):
        
    decodeChrom = ga.decodeChromosomes(lengths,chromosomes, boudList) #解码
    
    fitness,fitnessSum = ga.getFitness(ga.fitFunc1, decodeChrom)   #交叉变异选择前适应度
    
    newSelPopu = ga.SelNewPopu(chromosomes, fitnessSum) #轮赌选择新的群体
    
    crossoverPopu = ga.crossover(newSelPopu, pc = 0.8)  #单点交叉
    
    mutationPopu = ga.mutation(crossoverPopu, pm=0.003) #变异
    
    final_decode_chromo = ga.decodeChromosomes(lengths, mutationPopu, boudList)  #解码后的染色体
    
    decodeFitness,decodeFitnessSum = ga.getFitness(ga.fitFunc1, decodeChrom)    #进过交叉变异选择的适应度
    
    mergeChrom = np.concatenate((chromosomes, mutationPopu))   #合并两个原始群体和交叉变异之后的群体
    mergeFitness = np.concatenate((fitness, decodeFitness))   #合并适应度
    sortIndex = mergeFitness.argsort(axis=0)
    
    sortIndex  = [np.int(x) for x in sortIndex]
    sortIndex.reverse()
    
    chromosomes = mergeChrom[sortIndex[1:N], :]
#    chromosomes = mutationPopu.copy()
    
    maxFit = np.max(mergeFitness)
    optiVals.append(maxFit)
    index = np.where(mergeFitness == maxFit)
    optiSolutions.append(final_decode_chromo[index[0][0],:])
    print("iter：%d---适应度是：%f\n"%(it, maxFit) )


optiVal = np.max(optiVals)
index = np.where(optiVals == optiVal)
optiSolution = optiSolutions[index[0][0]]


plt.figure(1)
plt.plot(optiVals, c='r', label = 'Fitness')

plt.legend()
plt.show()

end = time.time()

print("目前找到的最优解是 x1:%f   x2:%f\n"%(optiSolution[0], optiSolution[1]))
print("最有解对应的适应度是：%f\n"%optiVal)
print("程序运行时间：%d s" %(end-start))








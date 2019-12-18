# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 09:08:37 2019

@author: Administrator
"""
import numpy as np
import matplotlib.pyplot as plt
import time
import NSGA_II as NSGA
import GA
import MOEAD

nowtim = time.time()
stim = time.time()
#-----------------------执行程序----------------------------
genration_num = 200  #迭代次数
n = 30                 #变量维数
Func_Num = 5          #进行测试的函数编号
Func = [NSGA.SCH, NSGA.FON, NSGA.KUR, MOEAD.ZDT1, MOEAD.ZDT2, MOEAD.ZDT3]     #进行测试的函数
bounds = [[-1000,1000], [-4,4], [-5,5], [0,1], [0,1], [0,1]] #测试函数对应的自变量边界
boudList = []         #变量边界
N = 100                #种群规模
delta = 0.01      #编码精度

#--------------------------------------------------------------------------------------
for i in range(n):
    boudList.append(bounds[Func_Num])  #初始化自变量边界


encodeLengths = GA.getEncodedLength(delta, boudList)   #获得编码的长度

Population = np.zeros((N , np.sum(encodeLengths) ))  
for i in range(N):
    Population[i,:] = np.random.randint(0,2,np.sum(encodeLengths))   #初始化种群 二进制编码
Population2 = Population.copy()     #初始化后代
vals = np.zeros((2*N,2))            #目标空间是2维
for iteration in range(genration_num):   #迭代次数
    mergePopu = np.concatenate( (Population, Population2) )     #父代子代群体混合
    decodePopu = GA.decodedChromosome(encodeLengths, mergePopu, boudList, delta)  #二进制->十进制
    vals = Func[Func_Num](decodePopu)
    
    t1 = time.time()
    front,rank = NSGA.fast_non_dominated_sort(vals)    #快速非支配排序
#    print('%f'%(time.time()-t1))
    P = []
    i = 0
    while (len(P)+len(front[i]))<=N:              #精英策略
        P.extend(front[i])
        if len(P)<N:
            i += 1
        else:
            break
    
    iDis  = NSGA.crowding_distance_assignment(vals, front[i])
    
    if(len(P)<N):
        index = np.argsort(-iDis, axis=0)
        #########################################
        valsIdx = np.argsort(vals[front[i], 0], axis=0)   #按第一个函数值进行排序  
        frontSort  = [front[i][valsIdx[k]] for k in range(len(valsIdx))]
        #######################################
        for k in range(N-len(P)):
             P.append(frontSort[index[k][0]])   #按照精英策略把剩下的部分个体添加  凑够N个个体
    Population = mergePopu[P]            #生成父代群体
    crossoverPopu = GA.crossover(Population, Pc=0.9)  #单点交叉
    mutationPopu = GA.mutation(crossoverPopu, Pm=1/np.sum(encodeLengths))  #位点变异
    Population2 = mutationPopu.copy() 
    
    if time.time() - nowtim>=3:
        nowtim = time.time()
        print("NSGA-II--%d/%d次--val（%f,%f）--mean（%f,%f）--%.3fs\n"
              %(iteration, genration_num, vals[0,0], vals[0,1],np.mean(vals[:N,0]),np.mean(vals[:N,1]),time.time()-stim))
        
    
decodePopuDisplay = GA.decodedChromosome(encodeLengths, Population, boudList, delta)  #二进制->十进制
valsDisplay =Func[Func_Num](decodePopuDisplay)


plt.scatter(valsDisplay[:,0], valsDisplay[:, 1], c='r', label='Parato Front')
plt.xlabel('func1')
plt.ylabel('func2')
plt.title('NSGA-II')
plt.legend()
plt.show()


#-------------------------保存数据------------------------
if   Func[Func_Num]==MOEAD.ZDT1:
    np.savetxt('ZDT1.txt', valsDisplay, fmt='%10.6f')
elif Func[Func_Num]==MOEAD.ZDT2:
    np.savetxt('ZDT2.txt', valsDisplay, fmt='%10.6f')
elif Func[Func_Num]==MOEAD.ZDT3:
    np.savetxt('ZDT3.txt', valsDisplay, fmt='%10.6f')



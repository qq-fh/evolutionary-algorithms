# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 13:01:08 2019

@author: Administrator
"""

import MOEAD
import matplotlib.pyplot as plt
import numpy as np

N = 70
T = 5
EP, epNum, RunningTime = MOEAD.MOEAD_Algorithm(MOEAD.ZDT3,   #测试函数
                                                   200,          #迭代次数
                                                   200,          #外部种群数量上限
                                                   N,            #内部种群数量
                                                   T,            #邻居数量
                                                   True          #是否打印运行数据
                                                   )    
x1 = np.array([x[0] for x in EP] ).reshape(len(EP),1)
x2 = np.array([x[1] for x in EP] ).reshape(len(EP),1)
vals = np.concatenate((x1,x2), axis=1 )
np.savetxt('ZDT3.txt', vals, fmt='%10.6f')

plt.scatter(x1, x2, color='r', label='PF' )
plt.xlabel('func1')
plt.ylabel('func2')
plt.title( 'MOEAD_ZDT1' )
plt.legend()          
plt.show()



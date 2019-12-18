# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 11:23:36 2019

@author: Administrator
"""

# !/usr/bin/env python
# -*- coding:utf-8 -*- 
import numpy as np
from scipy.optimize import fsolve, basinhopping
import random
import timeit
import matplotlib.pyplot as plt 

def getEncodedLength(delta=0.0001, boundarylist=[]):
	# 每个变量的编码长度
    lengths = []
    for i in boundarylist:
        lower = i[0]
        upper = i[1]
        # lamnda 代表匿名函数f(x)=0,50代表搜索的初始解
        res = fsolve(lambda x: ((upper - lower) * 1 / delta) - 2 ** x - 1, 50)
        length = int(np.floor(res[0]))
        lengths.append(length)
    return lengths


# 随机生成初始编码种群
def getIntialPopulation(encodelength, populationSize):
    # 随机化初始种群为0
    chromosomes = np.zeros((populationSize, sum(encodelength)), dtype=np.uint8)
    for i in range(populationSize):
        chromosomes[i, :] = np.random.randint(0, 2, sum(encodelength))
        # print('chromosomes shape:', chromosomes.shape)
    return chromosomes

# 染色体解码得到表现型的解
def decodedChromosome(encodelength, chromosomes, boundarylist, delta=0.0001):
    populations = chromosomes.shape[0]
    variables = len(encodelength)
    decodedvalues = np.zeros((populations, variables))
    for k, chromosome in enumerate(chromosomes):
        chromosome = chromosome.tolist()
        start = 0
        for index, length in enumerate(encodelength):
            # 将一个染色体进行拆分，得到染色体片段
            power = length - 1
            # 解码得到的10进制数字
            demical = 0
            for i in range(start, length + start):
                demical += chromosome[i] * (2 ** power)
                power -= 1
            lower = boundarylist[index][0]
            upper = boundarylist[index][1]
            decodedvalue = lower + demical * (upper - lower) / (2 ** length - 1)
            decodedvalues[k, index] = decodedvalue
            # 开始去下一段染色体的编码
            start += length
    return decodedvalues



# 得到个体的适应度值及每个个体被选择的累积概率
def getFitnessValue(func, chromosomesdecoded):
    # 得到种群规模和决策变量的个数
    population, nums = chromosomesdecoded.shape
    # 初始化种群的适应度值为0
    fitnessvalues = np.zeros((population, 1))
    # 计算适应度值
    for i in range(population):
        fitnessvalues[i, 0] = func(chromosomesdecoded[i, :])
        # 计算每个染色体被选择的概率
    probability = fitnessvalues / np.sum(fitnessvalues)
	# 得到每个染色体被选中的累积概率
    cum_probability = np.cumsum(probability)
    return fitnessvalues, cum_probability




# 新种群选择
def selectNewPopulation(chromosomes, cum_probability):
    m, n = chromosomes.shape
    newpopulation = np.zeros((m, n), dtype=np.uint8)
    # 随机产生M个概率值
    randoms = np.random.rand(m)
    for i, randoma in enumerate(randoms):
        logical = cum_probability >= randoma
        index = np.where(logical == 1)
        # index是tuple,tuple中元素是ndarray
        newpopulation[i, :] = chromosomes[index[0][0], :]
    return newpopulation


# 新种群交叉
def crossover(population, Pc=0.8):
    """
    :param population: 新种群
    :param Pc: 交叉概率默认是0.8
    :return: 交叉后得到的新种群
    """
    # 根据交叉概率计算需要进行交叉的个体个数
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
        # 随机产生一个交叉点
        crossoverPoint = random.sample(range(1, n), 1)
        crossoverPoint = crossoverPoint[0]
        # one-single-point crossover
        updatepopulation[a, 0:crossoverPoint] = population[a, 0:crossoverPoint]
        updatepopulation[a, crossoverPoint:] = population[b, crossoverPoint:]
        updatepopulation[b, 0:crossoverPoint] = population[b, 0:crossoverPoint]
        updatepopulation[b, crossoverPoint:] = population[a, crossoverPoint:]
    return updatepopulation



# 染色体变异
def mutation(population, Pm=0.01):
    """
    :param population: 经交叉后得到的种群
    :param Pm: 变异概率默认是0.01
    :return: 经变异操作后的新种群
    """
    updatepopulation = np.copy(population)
    m, n = population.shape
    # 计算需要变异的基因个数
    gene_num = np.uint8(m * n * Pm)
    # 将所有的基因按照序号进行10进制编码，则共有m*n个基因
    # 随机抽取gene_num个基因进行基本位变异
    mutationGeneIndex = random.sample(range(m * n), gene_num)
    # 确定每个将要变异的基因在整个染色体中的基因座(即基因的具体位置)
    for gene in mutationGeneIndex:
        # 确定变异基因位于第几个染色体
        chromosomeIndex = gene // n
        # 确定变异基因位于当前染色体的第几个基因位
        geneIndex = gene % n
        # mutation
        if updatepopulation[chromosomeIndex, geneIndex] == 0:
            updatepopulation[chromosomeIndex, geneIndex] = 1
        else:
            updatepopulation[chromosomeIndex, geneIndex] = 0
    return updatepopulation

# 定义适应度函数
def fitnessFunction():
    return lambda x: 21.5 + x[0] * np.sin(4 * np.pi * x[0]) + x[1] * np.sin(20 * np.pi * x[1])

def fitness_func1(x):
    n = len(x)
    S = 0
    for i in range(n):
        S += x[i]*np.sin(np.sqrt(abs(x[i])))
    return S
def fitness_func2(x):
    S = np.sum( x**2 - 10*np.cos(2*np.pi*x) + 10)
    return S
def fitness_func3(x):
    n = len(x)
    S = -20*np.exp(-0.2*np.sqrt(np.sum( x**2)/n) ) - np.exp( np.sum(np.cos(2*np.pi*x))/n ) + 20 + np.e
    return S

def main(max_iter=10000):
    #初始种群数量大小
    N = 10
    n = 30
    # 每次迭代得到的最优解
    optimalSolutions = []
    optimalValues = []
    # 决策变量的取值范围
    decisionVariables = []
    for i in range(n):
        decisionVariables.append([-5.12,5.12])
    # 得到染色体编码长度
    lengthEncode = getEncodedLength(boundarylist=decisionVariables)
    chromosomesEncoded = getIntialPopulation(lengthEncode, N)  #初始化种群编码
    for iteration in range(max_iter):
        # 种群解码
        decoded = decodedChromosome(lengthEncode, chromosomesEncoded, decisionVariables)
        # 得到个体适应度值和个体的累积概率
        evalvalues, cum_proba = getFitnessValue(fitness_func2, decoded)
        # 选择新的种群
        newpopulations = selectNewPopulation(chromosomesEncoded, cum_proba)
        # 进行交叉操作
        crossoverpopulation = crossover(newpopulations, Pc=0.7)
        # mutation
        mutationpopulation = mutation(crossoverpopulation, Pm=0.3)
        # 将变异后的种群解码，得到每轮迭代最终的种群
        final_decoded = decodedChromosome(lengthEncode, mutationpopulation, decisionVariables)
        # 适应度评价
        fitnessvalues, cum_individual_proba = getFitnessValue(fitness_func2, final_decoded)
        
        mergeChrom = np.concatenate((chromosomesEncoded, mutationpopulation))   #合并两个原始群体和交叉变异之后的群体
        mergeFitness = np.concatenate((evalvalues, fitnessvalues))   #合并适应度
#        statictic_decoded = decodedChromosome(lengthEncode, mergeChrom, decisionVariables)
        sortIndex = mergeFitness.argsort(axis=0)
        
        sortIndex  = [np.int(x) for x in sortIndex]
        
        chromosomesEncoded = mergeChrom[sortIndex[:N], :]
        
        
        # 搜索每次迭代的最优解，以及最优解对应的目标函数的取值
        minFit = np.min(list(mergeFitness))
        optimalValues.append(minFit)
#        index = np.where(mergeFitness == min(list(mergeFitness)))
#        optimalSolutions.append(statictic_decoded[index[0][0], :])
        print("iter：%d---fitness：%f"%(iteration, minFit))
        
        
    # 搜索最优解
    optimalValue = np.min(mergeFitness)
#    optimalIndex = np.where(mergeFitness == optimalValue)
#    optimalSolution = optimalSolutions[optimalIndex[0][0]]
    return optimalSolution, optimalValue,optimalValues
 
 
solution, value, optimalValues = main()

print('最优目标函数值:', value)
# 测量运行时间
plt.plot(optimalValues,  c='r', label='fitness')
plt.show()

 
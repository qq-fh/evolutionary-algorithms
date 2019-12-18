#ifndef _INC_H
#define _INC_H



#define N 300  //种群数量
#define T 20   //近邻数量
#define VAR_NUM 30 //输入变量维数
#define OBJ_NUM 2   //目标变量维数
#define K 500  //迭代次数

#define ZDT  0
#define DTLZ 1
#define F    2
#define FUNC_TYPE F
#define FUNC_NUM 2



#define VAR_MIN 0
#define VAR_MAX 1
#define LIMIT(x, min, max)  (x)>(max)?(max):((x)<(min)?(min):(x) )
#define pi 3.14159


void BubbleSort(double* data, int length, int* OutIndex );
void BubbleSort(double* data, int length);



#endif

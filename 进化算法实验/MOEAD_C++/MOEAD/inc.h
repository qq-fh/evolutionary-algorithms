#ifndef _INC_H
#define _INC_H



#define N 300  //��Ⱥ����
#define T 20   //��������
#define VAR_NUM 30 //�������ά��
#define OBJ_NUM 2   //Ŀ�����ά��
#define K 500  //��������

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

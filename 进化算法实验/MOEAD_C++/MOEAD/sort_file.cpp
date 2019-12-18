#include "sort_file.h"




//输入数组  返回升序下标
void BubbleSort(double* data, int length, int* OutIndex )
{
	for (int i=0;i<length;i++)
		OutIndex[i] = i;
	for (int i=0;i<length;i++)
	{
		for(int j=0;j<length-i-1;j++)
		{
			if(data[j]>data[j+1])
			{
				double temp = data[j];
				data[j] = data[j+1]; 
				data[j+1] = temp;
				int index_temp = OutIndex[j];
				OutIndex[j] = OutIndex[j+1];
				OutIndex[j+1] = index_temp;
			}
		}
	}
}
void BubbleSort(double* data, int length)
{
	for (int i=0;i<length;i++)
	{
		for(int j=0;j<length-i-1;j++)
		{
			if(data[j]>data[j+1])
			{
				double temp = data[j];
				data[j] = data[j+1]; 
				data[j+1] = temp;
			}
		}
	}
}
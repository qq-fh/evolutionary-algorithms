#include "iostream"
#include "MOEAD.h"
#include "ctime"
#include "random.h"
#include "fstream"
#include "time.h"


using namespace std;

double Z[OBJ_NUM] = {100,100};
double weight_vecs[N][OBJ_NUM];
int    neighbors[N][T];
double var_upper[VAR_NUM],var_lower[VAR_NUM];
Individual P[N];   //定义N个个体

void SavePF(double ObjTemp[][OBJ_NUM]);

int main()
{
	time_t start_t,end_t;
	start_t = clock();
	GenerateBoundary(var_lower, var_upper, 1);     //生成决策变量边界

	srand((int)time(NULL));

	cout << setiosflags(ios::fixed);
	cout.precision(3);

	int P_Idx[N];
	for (int i=0;i<N;i++)
	{
		P[i].SetNum(i);
		for(int j=0;j<VAR_NUM;j++)
			P[i].SetSolOnIndex(j, randomDouble(var_lower[j], var_upper[j],0.001));   //初始化解 Low-upper	精度0.001
		P_Idx[i] = i;
	}
	GetWeightVec();   //生成权重向量

	GetNeighbors();   //获取近邻下标信息

	double ObjVal[OBJ_NUM]={0};
	for(int i=0;i<N;i++)   //初始化Z
	{
		P[i].GetObjVal(ObjVal, FUNC_TYPE, FUNC_NUM);
		for(int j=0;j<OBJ_NUM;j++)
		{
			if(ObjVal[j] < Z[j])
				Z[j] = ObjVal[j];
		}
	}
	for (int Iter=0;Iter<K;Iter++)
	{
		double rd = randomDouble(0,1,0.01);
		if(rd>0.9)
			cout<<"Rate："<<(double)Iter/K*100.0<<"%"<<endl;
		for (int n=0; n<N; n++)  //主循环迭代
		{
			/***********************生成两个随机的下标***********************/
			double r = randomDouble(0, 1.0, 0.001);
			int Bound;
			if(r<DELTA)
				Bound = T;			//近邻中选取
			else
				Bound = N;          //整个种群中选取
			int Idx1 = randomInt(0, Bound);  
			int Idx2;
			while(1)
			{
				Idx2 = randomInt(0, Bound);   //生成两个不等的随机下标
				if (Idx1 != Idx2)
					break;
			}
			/*********************************DE**************************************/
			Individual newUnit;
			int *P_Para = new int[Bound]; //动态内存分配  
			if(Bound==T)
			{
				newUnit = P[n].DE(  P[ neighbors[n][Idx1] ], P[ neighbors[n][Idx2] ]  );    //DE
				memcpy(P_Para, neighbors[n], sizeof(int)*T);   //P为繁殖范围   选择是在近邻中、还是整个种群中
			}
			else if(Bound==N)
			{
				newUnit = P[n].DE(P[ Idx1 ], P[ Idx2 ]);
				memcpy(P_Para, P_Idx, sizeof(int)*N);    
			}
			else
				cout<<"Bound Err!!"<<endl;
			/**********************************更新Z**********************************/
			double newUnitObjVal[OBJ_NUM];
			newUnit.GetObjVal(newUnitObjVal, FUNC_TYPE, FUNC_NUM);  //获取对应的目标函数值
			for (int i=0;i<OBJ_NUM;i++)   //更新Z参考点   最小值
			{
				if(Z[i]>newUnitObjVal[i])
					Z[i] = newUnitObjVal[i];  
			}
			/***************************更新种群值（最多YITA个）********************************/
			int c=0;   //更新的值数量
			for (int i=0;i<Bound;i++)
			{
				if(c>=YITA)
					break;
				int Idx = randomInt(0, Bound-i);
				int J = P_Para[Idx];  //在P中随机选取一个数
				//切比雪夫聚合
				if(newUnit.GetTchebycheffVal(weight_vecs[J], Z) < P[J].GetTchebycheffVal(weight_vecs[J], Z)) 
				{
					for(int j=0;j<VAR_NUM;j++)
						P[J].SetSolOnIndex(j, newUnit.GetSolOnIndex(j) );    //更新J个体
					c++;
				}		
				for(int k=Idx;k<Bound-i;k++)
					P_Para[k] = P_Para[k+1];   //删除随机选到的解
			}
		}
	}
	double ObjTemp[N][OBJ_NUM];
	for (int i=0;i<N;i++)
	{
		P[i].GetObjVal(ObjTemp[i], FUNC_TYPE, FUNC_NUM);   //计算最终得到的解
		/*for(int j=0;j<OBJ_NUM;j++)
		cout<<ObjTemp[i][j]<<"  ";
		cout<<endl;*/
	}
	
	end_t = clock();
	SavePF(ObjTemp);    //保存PF
	/*for (int i=0;i<N;i++)
	P[i].PrintIndInfo();*/

	cout<<"程序运行总时间："<<(double)(end_t-start_t)/CLOCKS_PER_SEC<<"s"<<endl;
	  return 0;
}

void SavePF(double ObjTemp[][OBJ_NUM])
{
	ofstream outfile_pf;
	ofstream outfile_sol;
	if (FUNC_TYPE==ZDT)
	{
		switch(FUNC_NUM)
		{
		case 1:outfile_pf.open("ZDT1.txt");break;
		case 2:outfile_pf.open("ZDT2.txt");break;
		case 3:outfile_pf.open("ZDT3.txt");break;
		default:outfile_pf.open("default.txt");
		}
	}
	else if(FUNC_TYPE==DTLZ)
	{

	}
	else if(FUNC_TYPE==F)
	{
		switch(FUNC_NUM)
		{
		case 1:
			outfile_pf.open("F1.txt");
			outfile_sol.open("F1_SOL.txt");
			break;
		case 2:
			outfile_pf.open("F2.txt");
			outfile_sol.open("F2_SOL.txt");
			break;
		case 3:
			outfile_pf.open("F3.txt");
			outfile_sol.open("F3_SOL.txt");
			break;
		case 4:
			outfile_pf.open("F4.txt");
			outfile_sol.open("F4_SOL.txt");
			break;
		default:
			outfile_pf.open("default.txt");
		}
	}
	for (int i=0;i<N;i++)
	{
		outfile_pf<<ObjTemp[i][0]<<"  "<<ObjTemp[i][1];
		outfile_pf<<endl;
		for (int j=0;j<VAR_NUM;j++)
		{
			outfile_sol<<P[i].GetSolOnIndex(j)<<" ";
		}
		outfile_sol<<endl;
	}
	     
	outfile_pf.close();
	outfile_sol.close();
}
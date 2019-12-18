#ifndef _MOEAD_H
#define _MOEAD_H

#include "inc.h"
#include "iostream"
#include "iomanip"
#include "cmath"
#include "random.h"
using namespace  std;


extern double var_upper[VAR_NUM],var_lower[VAR_NUM];

class Individual   //个体类
{
private:
	int    m_num;            //序号
	double m_sol[VAR_NUM];   //个体对应的解
public:
	Individual();      //构造函数
	Individual(const Individual& obj);
	void SetNum(int num);
	int  GetNum(void);         //获得对应的序号
	void SetSol(double *sol, int var=VAR_NUM);   //设置个体对应的解
	void PrintIndInfo(void);
	void GetZdtVal(double* objVal, int num);
	void GetDtlzVal(double* objVal, int num);
	void GetF1_5Val(double* objVal, int num);
	void GetObjVal(double* objVal,int type, int num);   //获得目标空间函数值
	double GetSolOnIndex(int index);
	void SetSolOnIndex(int index, double val);
	bool isDominated(Individual obj);
	double GetTchebycheffVal(double *weightVec, double *Z);
	Individual DE(Individual p1, Individual p2);


	Individual operator*(const Individual& p)
	{
		Individual newUnit, newUnit1, newUnit2;
		////////////////SBX////////////////
		double yita1 = 1.0;
		double u = 0;
		double r = 0.0;
		double a,b;
		for (int i=0;i<VAR_NUM;i++)
		{
			u = randomDouble(0,1, 0.001);
			if (u<0.5)
				r = pow(2*u, 1.0/(yita1+1));
			else
				r = pow( 1/(2*(1-u)) , 1.0/(yita1+1));
			a = this->m_sol[i];
			b = p.m_sol[i];
		    newUnit1.m_sol[i] = 0.5*( (1+r)*a + (1-r)*b );
			newUnit2.m_sol[i] = 0.5*( (1-r)*a + (1+r)*b ); 
			newUnit1.m_sol[i]= LIMIT( newUnit1.m_sol[i], var_lower[i], var_upper[i]);    //自变量范围限制
			newUnit2.m_sol[i]= LIMIT( newUnit2.m_sol[i], var_lower[i], var_upper[i]);
		}
		if( newUnit1.isDominated(newUnit2) )
		{
			newUnit = newUnit2;
			//cout<<"Unit2 is Returned!"<<endl;
		}
		else
		{
			newUnit = newUnit1;
			//cout<<"Unit1 is Returned!"<<endl;
		}
		///////////////多项式变异////////////////
		double yita2 = 1.0;
		double delta = 0.0;
		for(int i=0;i<VAR_NUM;i++)
		{
			u = randomDouble(0, 1, 0.001);
			if(u>1.0/VAR_NUM)
				break;
			u = randomDouble(0, 1.0, 0.001);
			if(u<0.5)
				delta = pow(2*u, 1/(yita2+1.0))-1;
			else
				delta = 1 - pow( 2*(1-u), 1/(yita2+1.0) );
			newUnit.m_sol[i] += delta;
			newUnit.m_sol[i] = LIMIT(newUnit.m_sol[i], var_lower[i], var_upper[i]);
		}
		
		return newUnit;
	}

};
Individual::Individual()//构造函数
{
	m_num = 0;
	for(int i=0;i<VAR_NUM;i++)
		m_sol[i] = 0.0;
}   
Individual::Individual(const Individual& obj)   //拷贝构造函数
{
	memcpy(this->m_sol, obj.m_sol, sizeof(double)*VAR_NUM);
	this->m_num = obj.m_num;
}
double Individual::GetSolOnIndex(int index){return m_sol[index];}
void Individual::SetSolOnIndex(int index, double val){m_sol[index]=val;}
int  Individual::GetNum(void){return m_num;}
void Individual::SetNum(int num){m_num = num;}
void Individual::SetSol(double *sol, int var){memcpy(m_sol, sol, var*sizeof(double));}
Individual Individual::DE(Individual p2, Individual p3)   //差分算子
{
	Individual p1;
	/////////////////DE////////////////////
	double temp = 0.0;
	for (int i=0;i<VAR_NUM;i++)
	{
		temp = this->GetSolOnIndex(i)+F_PARA*( p2.GetSolOnIndex(i)-p3.GetSolOnIndex(i) );
		p1.SetSolOnIndex(i, LIMIT(temp, var_lower[i], var_upper[i] ) );
	}
	///////////////多项式变异////////////////
	double yita2 = 20.0;
	double delta = 0.0;
	double u=0.0;
	for(int i=0;i<VAR_NUM;i++)
	{
		u = randomDouble(0, 1, 0.001);
		if(u>1.0/VAR_NUM)
			continue;
		u = randomDouble(0, 1.0, 0.001);
		if(u<0.5)
			delta = pow(2*u, 1/(yita2+1.0))-1;
		else
			delta = 1 - pow( 2*(1-u), 1/(yita2+1.0) );
		double val = p1.GetSolOnIndex(i);
		p1.SetSolOnIndex(i, LIMIT(val+delta, var_lower[i], var_upper[i])  );
	}
	return p1;

}
void Individual::GetZdtVal(double* objVal, int num)
{
	double g = 0;
	for(int i=1;i<VAR_NUM;i++)
		g+=m_sol[i];
	g *= 9;
	g /= (VAR_NUM-1);
	g += 1;
	objVal[0] = m_sol[0];
	try
	{
		if (num==1)
			objVal[1] = g* (1- sqrt(objVal[0]/g));
		else if(num==2)
			objVal[1] = g* (1- pow(objVal[0]/g, 2.0));
		else if(num==3)
			objVal[1] = g* (1- sqrt(objVal[0]/g) - objVal[0]/g*sin(10*pi*m_sol[0]));
		if (g==0)
			throw"计算目标函数值时g=0";
	}catch(const char*msg)
	{
		cerr<<msg<<endl;
		objVal[1] = 0;
	}
}
void Individual::GetDtlzVal(double* objVal, int num)
{
	double g = 0.0;

}
void Individual::GetF1_5Val(double* objVal, int num)
{
	double f1=0.0;
	double f2=0.0;
	int nOdd = 0;
	int nEven = 0;
	for (int i=1;i<VAR_NUM;i++)
	{
		if( (i+1)/2==1 )//奇数
		{
			switch(num)
			{
				case 1:
					f1 += pow( m_sol[i] - pow(m_sol[0], 0.5*(1.0+3*(i-1)/(VAR_NUM-2))) , 2);
					break;
				case 2:
					f1 += pow( m_sol[i] - sin(6*pi*m_sol[0]+(i+1)*pi/VAR_NUM) , 2);
					break;
				case 3:
					f1 += pow( m_sol[i] - 0.8*m_sol[0]*cos(6*pi*m_sol[0]+(i+1)*pi/VAR_NUM) , 2);
					break;
				case  4:
					f1 += pow( m_sol[i] - 0.8*m_sol[0]*cos(  (6*pi*m_sol[0]+(i+1)*pi/VAR_NUM)/3.0) , 2);
					break;
			}
			nOdd++;
		}
		else		   //偶数
		{	
			switch(num)
			{
				case 1:
					f2 += pow( m_sol[i] - pow(m_sol[0], 0.5*(1.0+3*(i-1)/(VAR_NUM-2))) , 2);
					break;
				case 2:
					f2 += pow( m_sol[i] - sin(6*pi*m_sol[0]+(i+1)*pi/VAR_NUM) , 2);
					break;
				case 3:														
					f2 += pow( m_sol[i] - 0.8*m_sol[0]*sin(6*pi*m_sol[0]+(i+1)*pi/VAR_NUM) , 2);
					break;
				case 4:														
					f2 += pow( m_sol[i] - 0.8*m_sol[0]*sin(6*pi*m_sol[0]+(i+1)*pi/VAR_NUM) , 2);
					break;
			}

			nEven++;
		}
	}
	f1 /= nOdd;
	f1 *= 2;
	f1 += m_sol[0];
	objVal[0] = f1;

	f2 /= nEven;
	f2 *= 2;
	f2  = f2+1-sqrt(m_sol[0]);	
	objVal[1] = f2;

}
void Individual::GetObjVal(double* objVal,int type, int num)
{
	if (type==ZDT)
	{
		this->GetZdtVal(objVal, num);
	}
	else if(type==DTLZ)
	{
		this->GetDtlzVal(objVal, num);
	}
	else if(type==F)
	{
		this->GetF1_5Val(objVal, num);
	}
}
bool Individual::isDominated(Individual obj)
{
	double func_val1[OBJ_NUM],func_val2[OBJ_NUM];
	this->GetObjVal(func_val1, FUNC_TYPE, FUNC_NUM);
	obj.GetObjVal(func_val2, FUNC_TYPE, FUNC_NUM);
	bool noEqualFlag = false;
	bool existBiggerFlag = false;
	for(int i=0;i<OBJ_NUM;i++)
	{
		if( func_val1[i] < func_val2[i] )
			noEqualFlag = true;
		if(func_val1[i]>func_val2[i])
			existBiggerFlag = true;
	}
	if ( (noEqualFlag==true) || (existBiggerFlag==false) )
		return false;
	else 
		return true;
}
double Individual::GetTchebycheffVal(double *weightVec, double *Z)
{
	double vals[OBJ_NUM];
	double ObjVal[OBJ_NUM];
	this->GetObjVal(ObjVal, FUNC_TYPE, FUNC_NUM);
	for (int i=0;i<OBJ_NUM;i++)
		vals[i] = weightVec[i]*abs(ObjVal[i]-Z[i]);
	BubbleSort(vals, OBJ_NUM);
	return vals[OBJ_NUM-1];
}

void Individual::PrintIndInfo(void)
{
	cout << setiosflags(ios::fixed);
	cout.precision(3);
	cout<<"序号："<<m_num<<" \t";
	for(int i=0;i<VAR_NUM;i++)
		cout<<m_sol[i]<<" ";
	cout<<endl;
}

/**********************************子函数定义去***************************************/
extern double Z[OBJ_NUM] ;
extern double weight_vecs[N][OBJ_NUM];
extern int    neighbors[N][T];
extern double var_upper[VAR_NUM],var_lower[VAR_NUM];


void GenerateBoundary(double* lower, double* upper,int k=0, int ndim=VAR_NUM)
{
	for (int i=0;i<k;i++)
	{
		lower[i] = 0;
		upper[i] = 1;
	}
	for (int i=k;i<ndim;i++)
	{
		lower[i] = -1;
		upper[i] = 1;
	}
}

void GetWeightVec()
{
	double step = 1.0/N;
	for(int i=0;i<N;i++)
	{
		weight_vecs[i][0] = i*step;
		weight_vecs[i][1] = 1-i*step;
	}
}

void GetNeighbors()   //获得个体对应的近邻
{
	double dis[N][N];   //初始化为0
	int neighbor_index[N][N];
	for (int i=0;i<N;i++)
		for(int j=0;j<N;j++)
			dis[i][j] = 0.0;

	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++)
			for(int k=0;k<VAR_NUM;k++)
				dis[i][j] += pow(weight_vecs[i][k]-weight_vecs[j][k],2);
	for(int i=0;i<N;i++)
	{
		BubbleSort(dis[i],	N, neighbor_index[i]);
		memcpy(neighbors[i],neighbor_index[i], sizeof(int)*T ); //取出近邻下标
	}
}
#endif

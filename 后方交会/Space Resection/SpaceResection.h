#pragma once
#include"matrix.h"
#define _USE_MATH_DEFINES
#include<math.h>
#include<vector>
using namespace std;
#define PRECISION1 1.0e-3
#define PRECISION2 1.0e-6
struct portrayCoordinate//影像坐标
{
	double x,y;
};
struct floorCoordinate//地面坐标
{
	double X, Y,Z;
};
class SpaceResection
{
public://变量
	double m = 0,//像片比例尺
		f = 0,//焦距
		x0 = 0, y0 = 0;//内方位元素
	vector<portrayCoordinate>pic_XY;//像点坐标
	vector<floorCoordinate>floor_XY;//地面控制点坐标
	double pointnum=0;//地面控制点数量

	vector<double> M;// 保存六个值的中误差
	double	 m0 = 0,	// 单位权中误差
			vv = 0;	    // [vv]，即平方和

	double Xs = 0.0, Ys = 0.0, Zs = 0.0, t = 0.0, w = 0.0, k = 0.0;//外方位元素


public://函数
	SpaceResection();
	SpaceResection(double m, double f, double x0, double y0);
	
	void ReadCoordinate();//读取影像、地面坐标
	void GetStart();//初始化

	Matrix constructSR_R_Matrix(double a, double b, double c);//旋转矩阵R构造
	Matrix constructSR_A_Matrix(Matrix R, vector<double>&X, vector<double>&Y, vector<double>&Z);//矩阵A构造
	Matrix constructSR_L_Matrix(vector<double>X, vector<double>Y, vector<double>Z); // 矩阵L构造
	void correction(Matrix XX);//改正数
	bool CheckPrecison(Matrix &X);//收敛判断
	void AccuracyEvaluation(Matrix ATA, Matrix A, Matrix XX, Matrix L);//精度评定
	

	void calculate();//开始计算

	~SpaceResection();
};




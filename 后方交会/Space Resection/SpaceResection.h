#pragma once
#include"matrix.h"
#define _USE_MATH_DEFINES
#include<math.h>
#include<vector>
using namespace std;
#define PRECISION1 1.0e-3
#define PRECISION2 1.0e-6
struct portrayCoordinate//Ӱ������
{
	double x,y;
};
struct floorCoordinate//��������
{
	double X, Y,Z;
};
class SpaceResection
{
public://����
	double m = 0,//��Ƭ������
		f = 0,//����
		x0 = 0, y0 = 0;//�ڷ�λԪ��
	vector<portrayCoordinate>pic_XY;//�������
	vector<floorCoordinate>floor_XY;//������Ƶ�����
	double pointnum=0;//������Ƶ�����

	vector<double> M;// ��������ֵ�������
	double	 m0 = 0,	// ��λȨ�����
			vv = 0;	    // [vv]����ƽ����

	double Xs = 0.0, Ys = 0.0, Zs = 0.0, t = 0.0, w = 0.0, k = 0.0;//�ⷽλԪ��


public://����
	SpaceResection();
	SpaceResection(double m, double f, double x0, double y0);
	
	void ReadCoordinate();//��ȡӰ�񡢵�������
	void GetStart();//��ʼ��

	Matrix constructSR_R_Matrix(double a, double b, double c);//��ת����R����
	Matrix constructSR_A_Matrix(Matrix R, vector<double>&X, vector<double>&Y, vector<double>&Z);//����A����
	Matrix constructSR_L_Matrix(vector<double>X, vector<double>Y, vector<double>Z); // ����L����
	void correction(Matrix XX);//������
	bool CheckPrecison(Matrix &X);//�����ж�
	void AccuracyEvaluation(Matrix ATA, Matrix A, Matrix XX, Matrix L);//��������
	

	void calculate();//��ʼ����

	~SpaceResection();
};




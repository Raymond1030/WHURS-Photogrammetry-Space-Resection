#include "SpaceResection.h"

SpaceResection::SpaceResection()
{
}
SpaceResection::SpaceResection(double m, double f,double x0,double y0)
{
    this->m = m;
    this->f = f;
    this->x0 = x0;
    this->y0 = y0;
}






void SpaceResection::ReadCoordinate()
{
    FILE* fp = fopen("影像坐标.txt", "r");
    if (!fp)
    {
        cout << "读取失败";
        return;
    }
    fscanf_s(fp, "x\ty\t\n");

    while (!feof(fp))
    {
        portrayCoordinate p_XY;
        fscanf_s(fp, "%lf\t%lf\t\n", &p_XY.x, &p_XY.y);
        pic_XY.push_back(p_XY);
    }
    fclose(fp);
    this->pointnum = this->pic_XY.size();

    FILE* fs = fopen("地面坐标.txt", "r");
    if (!fs)
    {
        cout << "读取失败";
        return;
    }
    fscanf_s(fs, "X\tY\tZ\t\n");

    while (!feof(fs))
    {
        
        floorCoordinate fl_Cor;
        fscanf_s(fp, "%lf\t%lf\t%lf\t\n", &fl_Cor.X, &fl_Cor.Y,&fl_Cor.Z);
        floor_XY.push_back(fl_Cor);
    }
    fclose(fs);
    cout << "******数据读取成功******" << endl << endl;
    cout << "\t影像坐标" << "\t\t\t\t地面坐标" << endl;
    cout << "   x"<< "\t\t" << "   y";
    cout << "\t\t" << "   X" << "\t\t" << "   Y" << "\t\t" << "   Z" << endl;
    for (int i = 0; i < pointnum; i++)
    {
        
        cout  << pic_XY[i].x << "\t\t" << pic_XY[i].y;
        cout << "\t\t" << floor_XY[i].X << "\t\t" << floor_XY[i].Y << "\t\t" << floor_XY[i].Z<<endl;
    }

}

void SpaceResection::GetStart()//初始化，计算外方位元素的初值
{
    this->Zs = this->m * this->f;
    for (int i = 0; i <pointnum; i++)
    {
        this->pic_XY[i].x /= 1000;//单位换算 mm->m
        this->pic_XY[i].y /= 1000;
        this->Xs += this->floor_XY[i].X;
        this->Ys += this->floor_XY[i].Y;
    }
    this->Xs = this->Xs / pointnum;
    this->Ys = this->Ys / pointnum;
    t = w = k = 0;
}


Matrix SpaceResection::constructSR_R_Matrix(double a, double b, double c)
{
    Matrix R(3, 3);
    double sinA = sin(a), cosA = cos(a),
        sinB = sin(b), cosB = cos(b),
        sinC = sin(c), cosC = cos(c);
    R(0, 0) = cosA * cosC - sinA * sinB * sinC;
    R(0, 1) = -cosA * sinC - sinA * sinB * cosC;
    R(0, 2) = -sinA * cosB;
    R(1, 0) = cosB * sinC;
    R(1, 1) = cosB * cosC;
    R(1, 2) = -sinB;
    R(2, 0) = sinA * cosC + cosA * sinB * sinC;
    R(2, 1) = -sinA * sinC + cosA * sinB * cosC;
    R(2, 2) = cosA * cosB;
    return R;
}


Matrix SpaceResection::constructSR_A_Matrix(Matrix R, vector<double> &X, vector<double> &Y, vector<double> &Z)
{

    Matrix A(pointnum*2, 6);
    for (int i = 0; i < pointnum; i++)
    {
   
        Z[i]=(R(0, 2) * (floor_XY[i].X - Xs) + R(1, 2) * (floor_XY[i].Y - Ys) + R(2, 2) * (floor_XY[i].Z - Zs));
        //像主点近似坐标
        X[i]=(this->x0 - (this->f * (R(0, 0) * (floor_XY[i].X - Xs) + R(1, 0) * (floor_XY[i].Y - Ys) + R(2, 0) * (floor_XY[i].Z - Zs))) / Z[i]);
        Y[i]=(this->y0 - (this->f * (R(0, 1) * (floor_XY[i].X - Xs) + R(1, 1) * (floor_XY[i].Y - Ys) + R(2, 1) * (floor_XY[i].Z - Zs))) / Z[i]);

        //偏导数值
        A(i * 2, 0) = (R(0, 0) * f + R(0, 2) * (pic_XY[i].x - x0)) / Z[i];//a11=1/Z * (a1*f+a3(x-x0))
        A(i * 2, 1) = (R(1, 0) * f + R(1, 2) * (pic_XY[i].x - x0)) / Z[i];//a12=1/Z * (b1*f+b3(x-x0))
        A(i * 2, 2) = (R(2, 0) * f + R(2, 2) * (pic_XY[i].x - x0)) / Z[i];//a13=1/Z * (c1*f+c3(x-x0))

        A(i * 2 + 1, 0) = (R(0, 1) * f + R(0, 2) * (pic_XY[i].y - y0)) / Z[i];//a21=1/Z * (a2*f+a3(y-y0))
        A(i * 2 + 1, 1) = (R(1, 1) * f + R(1, 2) * (pic_XY[i].y - y0)) / Z[i];//a22=1/Z * (b2*f+b3(y-y0))
        A(i * 2 + 1, 2) = (R(2, 1) * f + R(2, 2) * (pic_XY[i].y - y0)) / Z[i];//a23=1/Z * (c2*f+c3(y-y0))
    
        A(i * 2, 3) = (pic_XY[i].y - y0) * sin(w) - ((pic_XY[i].x - x0) / f * ((pic_XY[i].x - x0) * cos(k) - (pic_XY[i].y - y0) * sin(k)) + f * cos(k)) * cos(w);//a14
        A(i * 2, 4) = -f * sin(k) -( (pic_XY[i].x - x0) / f * ((pic_XY[i].x - x0) * sin(k) + (pic_XY[i].y - y0) * cos(k)));//a15
        A(i * 2, 5) = pic_XY[i].y - y0;//a16
     
        A(i * 2 + 1, 3) = -(pic_XY[i].x - x0) * sin(w) - ((pic_XY[i].y - y0) / f * ((pic_XY[i].x - x0) * cos(k) - (pic_XY[i].y - y0) * sin(k)) - f * sin(k)) * cos(w);
        A(i * 2 + 1, 4) = -f * cos(k) - ((pic_XY[i].y - y0) / f * ((pic_XY[i].x - x0) * sin(k) + (pic_XY[i].y - y0) * cos(k)));
        A(i * 2 + 1, 5) = -(pic_XY[i].x - x0);

   
    }
    return A;
}

Matrix SpaceResection::constructSR_L_Matrix(vector<double>X,vector<double>Y,vector<double>Z)
{
    Matrix L(pointnum * 2, 1);
    for (int i = 0; i < pointnum; i++)
    {
        L(i * 2, 0) = pic_XY[i].x - X[i];
        L(i * 2+1, 0) = pic_XY[i].y - Y[i];
    }
    return L;
}

void SpaceResection::correction(Matrix XX)
{
    Xs = Xs + XX(0, 0);
    Ys = Ys + XX(1, 0);
    Zs = Zs + XX(2, 0);
    t = t + XX(3, 0);
    w = w + XX(4, 0);
    k = k + XX(5, 0);
}
bool SpaceResection::CheckPrecison(Matrix& X)
{
    bool Boolean;
    Boolean = { fabs(X(0,0)) < PRECISION1 &&fabs(X(1,0))<PRECISION1&& fabs(X(2,0)) < PRECISION1
        && fabs(X(3,0)) < PRECISION2 && fabs(X(4,0)) < PRECISION2 && fabs(X(5,0)) < PRECISION2 };

    return Boolean;
}

void SpaceResection::calculate()
{
    Matrix XX(6, 1);
    Matrix ATA,ATL;
    Matrix A,L;
    int Count=0;//迭代次数
    cout << "******开始迭代******" << endl;
    do{
        Count++;
        if (Count ==30) {
            cout << "迭代次数超限，可能不收敛" << endl;
            break;
        }
        Matrix R=constructSR_R_Matrix(this->t, this->w, this->k);//构造R矩阵

        vector<double>X(pointnum), Y(pointnum), Z(pointnum);      //像主点近似坐标
        A=constructSR_A_Matrix(R,X,Y,Z);//构造A矩阵
        L=constructSR_L_Matrix(X,Y,Z);//构造L矩阵

        //X=inv(A^T *A) * A^T * L
        ATA = A.transpose() * A;
        ATL = A.transpose() * L;
        XX= ATA.inverse()* ATL;
        correction(XX);

        cout << "迭代次数第" << Count << "次" << endl;
        cout << "外方位元素" << endl;
        cout << "\tXs =" << Xs << "\tYs =" << Ys << "\tZs =" << Zs << endl;
        cout << "\tφ =" << t << "\tω = " << w << "\tκ =" << k << endl<<endl;
    } while (!CheckPrecison(XX));

    AccuracyEvaluation(ATA, A, XX, L);
    cout << "\n单位权中误差:" << m0<<endl;
    cout << "外方位元素的精度评定：" << endl;
    cout << "Xs:" << M[0] << endl;
    cout << "Ys:" << M[1] << endl;
    cout << "Zs:" << M[2] << endl;  
    cout << "φ:" << M[3] << endl;
    cout << "ω:" << M[4] << endl;
    cout << "κ:" << M[5] << endl;
}

void SpaceResection::AccuracyEvaluation(Matrix ATA,Matrix A,Matrix XX,Matrix L)
{
    // 精度评定
    vector<vector<double>> Q(6, vector<double>(1));
    for (int i = 0; i < 6; i++) {
        Q[i][0] = ATA(i, i);
    }


    Matrix V = A * XX - L;			// 当有 N 个控制点时：V = A × X - L

    for (int i = 0; i < 8; i++) {
        vv = vv + V(i, 0) * V(i, 0);
    }
    m0 = sqrt(vv / (2 * pointnum - 6));			// 单位全中误差 m0

    for (int i = 0; i < 6; i++) {
        double Qi = Q[i][0];
        M.push_back(m0 * sqrt(Qi));
        if (i > 2) {
            M[i] = M[i] * 180 * 3600 / M_PI;	// 转换为角度制
        }
    }
}

SpaceResection::~SpaceResection()
{
}
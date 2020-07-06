#pragma once
#include"GlobalData.h"
class Nfunction
{
protected:
	GlobalData d;
public:
	double* dNdksi;
	double* dNdeta;
	double* N;
	Nfunction();
	void setAllNfun(double ksi , double eta);
	void setAlldNfun(double ksi, double eta);
	double calcN1(double ksi, double eta);
	double calcN2(double ksi, double eta);
	double calcN3(double ksi, double eta);
	double calcN4(double ksi, double eta);
};


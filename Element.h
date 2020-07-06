#pragma once
#include"Node.h"
class Element {
protected:
	
	double ID[4];
public:
	Element();
	double getID0();
	double getID1();
	double getID2();
	double getID3();
	void SETElement(double, double, double, double);
};
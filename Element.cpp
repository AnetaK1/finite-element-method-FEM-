#include"Element.h"
#include<iostream>

using std::cout;
using std::endl;

Element::Element() {
	ID[0] = 0;
	ID[1] = 0;
	ID[2] = 0;
	ID[3] = 0;
}

double Element::getID0() {
	return ID[0];
}

double Element::getID1() {
	return ID[1];
}
double Element::getID2() {
	return ID[2];
}
double Element::getID3() {
	return ID[3];
}

//funkcja do ustwiania Id wêz³ów w elemencie
void Element::SETElement(double a1, double a2, double a3, double a4) {
	 
	ID[0]=a1 ;
	ID[1]=a2 ;
	ID[2]=a3 ;
	ID[3]=a4 ;
}


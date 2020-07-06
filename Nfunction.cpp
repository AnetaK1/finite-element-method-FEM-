#include "Nfunction.h"


Nfunction::Nfunction() {
	
	N = new double[d.get_NElemNodes()];
	dNdeta= new double[d.get_NElemNodes()];
	dNdksi= new double[d.get_NElemNodes()];

	for (int i = 0; i < d.get_NElemNodes(); i++) {
		N[i] = 0;
		dNdeta[i] = 0;
		dNdksi[i] = 0;
	}
}

double Nfunction::calcN1(double ksi, double eta) {
	return (0.25 * (1 - ksi) * (1 - eta));
}

double Nfunction::calcN2(double ksi, double eta) {
	return (0.25 * (1 + ksi) * (1 - eta));
}

double Nfunction::calcN3(double ksi, double eta) {
	return (0.25 * (1 + ksi) * (1 + eta));
}

double Nfunction::calcN4(double ksi, double eta) {
	return (0.25 * (1 - ksi) * (1 + eta));
}


void Nfunction::setAllNfun(double ksi,double eta) {
	N[0] = calcN1(ksi, eta);
	N[1] =calcN2(ksi,eta);
	N[2] = calcN3(ksi,eta);
	N[3] = calcN4(ksi,eta);
}


void Nfunction::setAlldNfun(double ksi, double eta) {
	dNdksi[0] = -(0.25 * (1 - eta));
	dNdksi[1] = (0.25 * (1 - eta));
	dNdksi[2] = (0.25 * (1 + eta));
	dNdksi[3] = -(0.25 * (1 + eta));
	dNdeta[0] = -(0.25 * (1 - ksi));
	dNdeta[1] = -(0.25 * (1 + ksi));
	dNdeta[2] = (0.25 * (1 + ksi));
	dNdeta[3] = (0.25 * (1 - ksi));

}
#include "Gauss_points.h"
#include<cmath>

//wprowadzenie wartoœci punktów ca³kowania
Gauss_points::Gauss_points() {
	war[0] = (1 / (sqrt(3)));
	war[1] = -(1 / (sqrt(3)));
	war[2] = -1;
	war[3] = 1;
	//punkty ca³kowania do ca³kowania po objêtoœci
	pc[0][0] = war[1];//punkt1ksi p[0][0]-drugie 0 to ksi
	pc[0][1] = war[1];//punkt1eta p[0][1]-1 to eta
	pc[1][0] = war[0];
	pc[1][1] = war[1];
	pc[2][0] = war[0];
	pc[2][1] = war[0];
	pc[3][0] = war[1];
	pc[3][1] = war[0];
	//punkty ca³kowania do ca³kowania po powierzchni
	pc_pow[0][0] = war[1];
	pc_pow[0][1] = war[2];
	pc_pow[1][0] = war[0];
	pc_pow[1][1] = war[2];
	pc_pow[2][0] = war[3];
	pc_pow[2][1] = war[1];
	pc_pow[3][0] = war[3];
	pc_pow[3][1] = war[0];
	pc_pow[4][0] = war[0];
	pc_pow[4][1] = war[3];
	pc_pow[5][0] = war[1];
	pc_pow[5][1] = war[3];
	pc_pow[6][0] = war[2];
	pc_pow[6][1] = war[0];
	pc_pow[7][0] = war[2];
	pc_pow[7][1] = war[1];
}

#pragma once
class GlobalData {	
protected:	
	double W, H, nW, nH, k, ro, c,alfa, t_amb,t_init,sim_time,sim_step;
	int nE, nN;
	int nElemNodes;

public:
	GlobalData();
	double get_W();
	double get_H();
	double get_nW();
	double get_nH();
	int get_nN();
	int get_nE();
	int get_NElemNodes();
	double get_k();
	double get_c();
	double get_ro();
	double get_alfa();
	double get_t_amb();
	double get_t_init();
	double get_sim_time();
	double get_sim_step();

	
};
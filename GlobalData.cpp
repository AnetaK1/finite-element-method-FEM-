#include"GlobalData.h"
#include<iostream>
#include<fstream>
#include<cstdlib>

using std::cout;
using std::endl;
using std::fstream;
using namespace std;

GlobalData::GlobalData() {
	
	//wczytywanie danych z pliku
	ifstream file;
	file.open("test_case1.txt");
	file >> t_init;
	file >> sim_time;
	file >> sim_step;
	file >> t_amb;
	file >> alfa;
	file >> H;
	file >> W;
	file >> nH;
	file >> nW;
	file >> c;
	file >> k;
	file >> ro;
	file >> nElemNodes;
	nN = ((nW) * (nH));
	nE = (nH - 1) * (nW - 1);
	
}


double GlobalData::get_W() {
	
	return this->W;
}

double GlobalData::get_H() {
	 
	return this->H;
}

double GlobalData::get_nW() {
	
	return this->nW;
}

double GlobalData::get_nH() {
	
	return  this->nH;
}

int GlobalData::get_nN() {
	return this->nN;
}

int GlobalData::get_nE() {
	return this->nE;
}

int GlobalData::get_NElemNodes() {
	return this->nElemNodes;
}

double GlobalData::get_k() {
	return this->k;
}

double GlobalData::get_ro() {
	return this->ro;
}

double GlobalData::get_c() {
	return this->c;
}

double GlobalData::get_alfa() {
	return this->alfa;
}

double GlobalData::get_t_amb() {
	return t_amb;
}

double GlobalData::get_sim_step() {
	return sim_step;
}

double GlobalData::get_sim_time() {
	return sim_time;
}

double GlobalData::get_t_init() {
	return t_init;
}
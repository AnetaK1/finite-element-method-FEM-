#pragma once
#include"Element.h"
#include"Node.h"
#include"GlobalData.h"
#include"Universal_Element.h"

class Grid {
protected:
	double dx, dy;
	Node* nd;
	Element* elem;
	Universal_Element* l;
	GlobalData d;
	double** H_matrix_agr;
	double** C_matrix_agr;
	double* P_vect_agr;
	double* t0_vect;
	double* t1_vect;
	void print_min_max_temp();
	
	double** mat;
	double* gaussianElimination(double** mat);
	void swap_row(double** mat, int i, int j);
	int forwardElim(double** mat);
	double* backSub(double** mat);


public:
	Grid();
	void create_agregation();
	void equations();
	void print_grid();
	void print_matrix();
	~Grid();
};
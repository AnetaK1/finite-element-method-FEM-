#pragma once
#include"GlobalData.h"
#include"Gauss_points.h"
#include"Nfunction.h"
#include"Element.h"
#include"Node.h"
class Universal_Element
{
protected:
	GlobalData d;


	double determinant;
	double det_jakobi_area;
	double jakobi[2][2];
	double jakobi_reverse[2][2];
	double** deria_dn_deta;
	double** deria_dn_dksi;
	double dndx[4][4];
	double dndy[4][4];
	double Nf[4][4];
	int* id;
	Nfunction* nFun;
	Gauss_points* pGauss;


	void create_jakobi(Element el, Node* nd, int k);
	void create_reverse_jakobi();
	double jakobi_det_area( Node nd1, Node nd2);
public:
	double C[4][4];
	double H[4][4];
	double Hbc[4][4];
	double P[4];
	Universal_Element();
	void print();

	void create_matrixes(Element el, Node* nd);
	void create_vector_P(Element el, Node* nd);
	~Universal_Element();

};


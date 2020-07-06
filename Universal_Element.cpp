#include "Universal_Element.h"
#include<cmath>
#include<iostream>

using namespace std;


Universal_Element::Universal_Element() {
	determinant = 0;

	id = new int[d.get_NElemNodes()];
	nFun = new Nfunction[d.get_NElemNodes()];
	pGauss = new Gauss_points();

	deria_dn_deta = new double* [d.get_NElemNodes()];
	deria_dn_dksi = new double* [d.get_NElemNodes()];
	for (int i = 0; i < d.get_NElemNodes(); i++) {
		deria_dn_deta[i] = new double[d.get_NElemNodes()];
		deria_dn_dksi[i] = new double[d.get_NElemNodes()];

	}

	

	//obliczenie wartosci funkcji kszta³tu oraz ich pochodnych po ksi oraz eta dla punktów ca³kowania po objêtoœci
	for (int i = 0; i < 4; i++) {
		nFun[i].setAlldNfun(pGauss->pc[i][0], pGauss->pc[i][1]);
		nFun[i].setAllNfun(pGauss->pc[i][0], pGauss->pc[i][1]);

	}


	//wpisanie wartoœci pochodnych funkcji kszta³tu po ksi oraz eta z wszystkich punktów ca³kowania do macierzy 
	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < 4; i++) {

			deria_dn_deta[i][j] = nFun[i].dNdeta[j];
			deria_dn_dksi[i][j] = nFun[i].dNdksi[j];
		}
	}



	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			jakobi[i][j] = 0;
			jakobi_reverse[i][j] = 0;
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Hbc[i][j] = 0;
			H[i][j] = 0;
			C[i][j] = 0;
		}
		P[i] = 0;
	}




}


void Universal_Element::create_jakobi(Element el, Node* nd, int k) {
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			jakobi[i][j] = 0;
			jakobi_reverse[i][j] = 0;
		}
	}
	id[0] = el.getID0();
	id[1] = el.getID1();
	id[2] = el.getID2();
	id[3] = el.getID3();
	
	//obliczanie macierzy jakobiego dla poszczególnych punktów ca³kowania, k- oznacza numer punktu ca³kowania
	for (int i = 0; i < 4; i++) {
		jakobi[0][0] += deria_dn_dksi[k][i] * (nd[id[i]].get_x());
		jakobi[0][1] += deria_dn_deta[k][i] * (nd[id[i]].get_x());
		jakobi[1][0] += deria_dn_dksi[k][i] * (nd[id[i]].get_y());
		jakobi[1][1] += deria_dn_deta[k][i] * (nd[id[i]].get_y());

	}
	//obliczenie wyznacznika 
	determinant = jakobi[0][0] * jakobi[1][1] - jakobi[0][1] * jakobi[1][0];

}


//tworzenie odwrotnej macierzy jakobiego
void Universal_Element::create_reverse_jakobi() {
	
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			jakobi_reverse[i][j] = jakobi[i][j] * (1 / determinant);
		}
	}
}



void Universal_Element::create_matrixes(Element el, Node* nd) {

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Hbc[i][j] = 0;
			H[i][j] = 0;
			C[i][j] = 0;
		}
		P[i] = 0;
	}
	//stworzenie wektora P oraz macierzy Hbc
	create_vector_P(el,nd);


	//pêtla po punktach ca³kownia
	for (int l = 0; l < 4; l++) {
		//tworzenie macierzy jakobiego, odwrotnej oraz wyznacznika dla poszcególnych punktów ca³kownaia
		create_jakobi(el, nd,l);
		create_reverse_jakobi();
		 
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j <4; j++) {

				//obliczanie pochodnych funkcji kszta³tu po x oraz y dla poszcególnych punktów ca³kowania
				dndx[l][j] = jakobi_reverse[0][0] * deria_dn_dksi[l][j] + jakobi_reverse[0][1] * deria_dn_deta[l][j];
				dndy[l][j] = jakobi_reverse[1][0] * deria_dn_dksi[l][j] + jakobi_reverse[1][1] * deria_dn_deta[l][j];
				//wpisanie wartoœci funkcji kszta³tu do macierzy przechowuj¹cej wartoœci funkcji kszta³tu dla ka¿dego punktu ca³kownia
				Nf[l][j] = nFun[l].N[j];
				//ca³kowanie macierzy H oraz C dla poszcególnych punktów ca³kowania i od razu je sumuje ze sob¹
				H[i][j] += ((dndx[l][i] * dndx[l][j] + dndy[l][i] * dndy[l][j]) * d.get_k() * determinant);
				C[i][j] += (Nf[l][i] * Nf[l][j]) * d.get_c() * d.get_ro() * determinant ;
				
			}

		}
	}

	//dodanie macierzy Hbc do mcierzy H
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[i][j] += Hbc[i][j];
		}
	}


}

//obliczenie jakobianu przekszta³cenia dla przypadku 1D, det[J]=L/2
double Universal_Element::jakobi_det_area( Node nd1,Node nd2) {
	double wynik= sqrt((pow(nd1.get_x() - nd2.get_x(),2)) + (pow(nd1.get_y() - nd2.get_y(),2)));
	wynik = wynik / 2;
	return  wynik;
	
}


void Universal_Element::create_vector_P(Element elem, Node* nd) {
	
	int id[4] = { elem.getID0(),elem.getID1(),elem.getID2(),elem.getID3() };


	double det = 0;
	double pc1[4];
	double pc2[4];
	// sprawdzanie czy dana œciana ma warunek brzegowy, jesli ma obliczam macierz Hbc oraz wektor P
	//obliczam je dla 2 punktów ca³kowania, które znajduj¹ siê na danej œcianie, sumuje 
	//w macierzy Hbc mno¿e wartoœci w punktach ca³kowania mno¿e jeszce przez alfa, wartosci dodane z dwóch punktów ca³kowania prze det[J]
	// w wektorze P robie podobnie jak w macierzy Hbc, jedynie w punktach ca³kowania domna¿am jeszcze prze temperature otoczenia
	if ((nd[id[0]].get_BC() == 1) && (nd[id[1]].get_BC() == 1)) {
		
		det = jakobi_det_area(nd[id[0]], nd[id[1]]);
		pc1[0]= nFun->calcN1(pGauss->pc_pow[0][0], pGauss->pc_pow[0][1]);
		pc1[1]= nFun->calcN2(pGauss->pc_pow[0][0], pGauss->pc_pow[0][1]);
		pc2[0]= nFun->calcN1(pGauss->pc_pow[1][0], pGauss->pc_pow[1][1]);
		pc2[1]= nFun->calcN2(pGauss->pc_pow[1][0], pGauss->pc_pow[1][1]);
		pc1[2] = pc1[3] = pc2[2] = pc2[3] = 0;
		
		 P[0] += -((d.get_alfa() * pc1[0]*d.get_t_amb()) + (d.get_alfa() * pc2[0]* d.get_t_amb()))*det;
		 P[1]+= -((d.get_alfa() * pc1[1]* d.get_t_amb()) + (d.get_alfa() * pc2[1] * d.get_t_amb()))*det;

		 Hbc[0][0]+=((pc1[0]*pc1[0]*d.get_alfa())+(pc2[0]*pc2[0]* d.get_alfa()))*det;
		 Hbc[0][1]+= ((pc1[0] * pc1[1] * d.get_alfa()) + (pc2[0] * pc2[1] * d.get_alfa())) * det;
		 Hbc[1][0]+= ((pc1[1] * pc1[0] * d.get_alfa()) + (pc2[1] * pc2[0] * d.get_alfa())) * det;
		 Hbc[1][1]+= ((pc1[1] * pc1[1] * d.get_alfa()) + (pc2[1] * pc2[1] * d.get_alfa())) * det;
		
	}
	 if ((nd[id[1]].get_BC() == 1) && (nd[id[2]].get_BC() == 1)) {
		det = jakobi_det_area(nd[id[1]], nd[id[2]]);
		pc1[0] = pc1[3] = pc2[0] = pc2[3] = 0;
		pc1[1]= nFun->calcN2(pGauss->pc_pow[2][0], pGauss->pc_pow[2][1]);
		pc1[2]= nFun->calcN3(pGauss->pc_pow[2][0], pGauss->pc_pow[2][1]);
		pc2[1]= nFun->calcN2(pGauss->pc_pow[3][0], pGauss->pc_pow[3][1]);
		pc2[2]= nFun->calcN3(pGauss->pc_pow[3][0], pGauss->pc_pow[3][1]);
		

		
		P[1]+= -((d.get_alfa() * pc1[1] * d.get_t_amb()) + (d.get_alfa() * pc2[1] * d.get_t_amb()))*det;
		P[2]+= -((d.get_alfa() * pc1[2] * d.get_t_amb()) + (d.get_alfa() * pc2[2] * d.get_t_amb()))*det;

		Hbc[1][1] += ((pc1[1] * pc1[1] * d.get_alfa()) + (pc2[1] * pc2[1] * d.get_alfa())) * det;
		Hbc[1][2] += ((pc1[1] * pc1[2] * d.get_alfa()) + (pc2[1] * pc2[2] * d.get_alfa())) * det;
		Hbc[2][1] += ((pc1[2] * pc1[1] * d.get_alfa()) + (pc2[2] * pc2[1] * d.get_alfa())) * det;
		Hbc[2][2] += ((pc1[2] * pc1[2] * d.get_alfa()) + (pc2[2] * pc2[2] * d.get_alfa())) * det;
		

	}
	if ((nd[id[2]].get_BC() == 1) && (nd[id[3]].get_BC() == 1)) {
		det = jakobi_det_area(nd[id[2]], nd[id[3]]);

		pc1[0] = pc1[1] = pc2[0] = pc2[1] = 0;
		pc1[2]= nFun->calcN3(pGauss->pc_pow[4][0], pGauss->pc_pow[4][1]);
		pc1[3]= nFun->calcN4(pGauss->pc_pow[4][0], pGauss->pc_pow[4][1]);
		pc2[2]= nFun->calcN3(pGauss->pc_pow[5][0], pGauss->pc_pow[5][1]);
		pc2[3]= nFun->calcN4(pGauss->pc_pow[5][0], pGauss->pc_pow[5][1]);
		
		
		P[2]+= -((d.get_alfa() * pc1[2] * d.get_t_amb()) + (d.get_alfa() * pc2[2] * d.get_t_amb()))*det;
		P[3]+= -((d.get_alfa() * pc1[3] * d.get_t_amb()) + (d.get_alfa() * pc2[3] * d.get_t_amb()))*det;

		Hbc[2][2] += ((pc1[2] * pc1[2] * d.get_alfa()) + (pc2[2] * pc2[2] * d.get_alfa())) * det;
		Hbc[2][3] += ((pc1[2] * pc1[3] * d.get_alfa()) + (pc2[2] * pc2[3] * d.get_alfa())) * det;
		Hbc[3][2] += ((pc1[3] * pc1[2] * d.get_alfa()) + (pc2[3] * pc2[2] * d.get_alfa())) * det;
		Hbc[3][3] += ((pc1[3] * pc1[3] * d.get_alfa()) + (pc2[3] * pc2[3] * d.get_alfa())) * det;
		
	}
	 if ((nd[id[3]].get_BC() == 1) && (nd[id[0]].get_BC() == 1)) {
		 
		det = jakobi_det_area(nd[id[3]], nd[id[0]]);
	
		pc1[1] = pc1[2] = pc2[1] = pc2[2] = 0;
		pc1[0]= nFun->calcN1(pGauss->pc_pow[6][0], pGauss->pc_pow[6][1]);
		pc1[3]= nFun->calcN4(pGauss->pc_pow[6][0], pGauss->pc_pow[6][1]);
		pc2[0]= nFun->calcN1(pGauss->pc_pow[7][0], pGauss->pc_pow[7][1]);
		pc2[3]= nFun->calcN4(pGauss->pc_pow[7][0], pGauss->pc_pow[7][1]);
		
		
		P[0]+= -((d.get_alfa() * pc1[0] * d.get_t_amb())+(d.get_alfa() *pc2[0] * d.get_t_amb()))*det;
		P[3]+= -((d.get_alfa() * pc1[3] * d.get_t_amb()) + (d.get_alfa() * pc2[3] * d.get_t_amb()))*det;

		Hbc[0][0] += ((pc1[0] * pc1[0] * d.get_alfa()) + (pc2[0] * pc2[0] * d.get_alfa())) * det;
		Hbc[0][3] += ((pc1[0] * pc1[3] * d.get_alfa()) + (pc2[0] * pc2[3] * d.get_alfa())) * det;
		Hbc[3][0] += ((pc1[3] * pc1[0] * d.get_alfa()) + (pc2[3] * pc2[0] * d.get_alfa())) * det;
		Hbc[3][3] += ((pc1[3] * pc1[3] * d.get_alfa()) + (pc2[3] * pc2[3] * d.get_alfa())) * det;
		
	
	}
	

}

void Universal_Element::print() {

	cout << "Wyznacznik: " << determinant << endl;

	cout << "dn_deta" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << deria_dn_deta[i][j] << " ";
		}
		cout << endl;
	}
	cout << "dn_dksi" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << deria_dn_dksi[i][j] << " ";
		}
		cout << endl;
	}



	cout << "Macierz jakobi:" << endl;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			cout << jakobi[i][j] << " ";
		}
		cout << endl;
	}

	cout << "Macierz jakobi odwrotna: " << endl;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			cout << jakobi_reverse[i][j] << " ";
		}
		cout << endl;
	}

	cout << "macierzH " << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << H[i][j] << "  ";
		}
		cout << endl;
	}

	cout << "macierzC " << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << C[i][j] << "  ";
		}
		cout << endl;
	}
	
	
	cout << "Macierz Bc" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << Hbc[i][j] << " ";
		}
		cout << endl;
	}

	cout << "Wektor P" << endl;
	for (int i = 0; i < 4; i++) {
		cout << P[i] << " ";
	}
	cout << endl;
	
}




Universal_Element::~Universal_Element() {
	for (int i = 0; i < d.get_NElemNodes(); i++) {
		delete[] deria_dn_deta[i];
		delete[] deria_dn_dksi[i];
	}
	delete[]deria_dn_deta;
	delete[]deria_dn_dksi;
	delete[]id;
	delete[] nFun;
	delete pGauss;
}

#include"Grid.h"
#include<iostream>
#include<fstream>

using namespace std;
using std::cout;
using std::endl;

Grid::Grid() {

	dx = d.get_W() / (d.get_nW() - 1);
	dy = d.get_H() / (d.get_nH() - 1);
	
	nd = new Node[d.get_nN()];
	elem = new Element[d.get_nE()];
	l = new Universal_Element[d.get_nE()];
	P_vect_agr = new double[d.get_nN()];
	H_matrix_agr = new double* [d.get_nN()];
	C_matrix_agr = new double* [d.get_nN()];
	t0_vect = new double[d.get_nN()];
	t1_vect = new double[d.get_nN()];
	mat = new double*[d.get_nN()];
	for (int i = 0; i < d.get_nN(); i++) {
		H_matrix_agr[i] = new double[d.get_nN()];
		C_matrix_agr[i] = new double[d.get_nN()];
		mat[i] = new double[d.get_nN() + 1];
	}

	for (int i = 0; i < d.get_nN(); i++) {
		for (int j = 0; j < d.get_nN(); j++) {
			H_matrix_agr[i][j] = 0;
			C_matrix_agr[i][j] = 0;

		}
		
		P_vect_agr[i] = 0;
	}

	


	//przypisywanie wartoœci wspó³rzêdnych x i y oraz tempertaury pocz¹tkowej oraz warunku brzegowego 
	//poszczególnym wêz³om
	int k = 0;
		for (int i = 0; i < d.get_nW(); i++) {
			for (int j = 0; j < d.get_nH(); j++) {
			
					nd[k].set_x(i * dx);
					nd[k].set_y( j * dy);
					nd[k].set_t(d.get_t_init());
					

					if ((nd[k].get_x() == 0 || nd[k].get_y() == 0) || (nd[k].get_x() == (d.get_nW() - 1 )*dx || nd[k].get_y() == (d.get_nH() - 1 )*dy)) {
						nd[k].set_BC(true) ;
					
					}
					
				k++;
				
			}
		
	    }
	
		
		//przypisanie odpowiednich id wêz³ów elementom
		int l = 0;
		for (int i = 0; i < d.get_nW()-1; i++) {
			for (int j = 0; j < d.get_nH()-1; j++) {
				      
					l = i * (d.get_nW() - 1) + j;
				elem[l].SETElement(l+i, (l+i + d.get_nH()), (l+i + d.get_nH() + 1), (l+i + 1));
			}

		}
	
	
}

void Grid::create_agregation() {
	
	//pêtla po wszystkich elemnetach w której tworzê lokalne macierze H, C oraz wektor P 
	//kolejnie wartoœci s¹ agregowane
	for (int k = 0; k < d.get_nE(); k++) {
		int id[4] = { elem[k].getID0(),elem[k].getID1(),elem[k].getID2(),elem[k].getID3() };

		l[k].create_matrixes(elem[k], nd);
		
		for (int i = 0; i <d.get_NElemNodes(); i++) {

				P_vect_agr[id[i]] += l[k].P[i];
			
			for (int j = 0; j < d.get_NElemNodes(); j++) {
			
					H_matrix_agr[id[i]][id[j]] += l[k].H[i][j];
					C_matrix_agr[id[i]][id[j]] += l[k].C[i][j];
				
			}
			
			
		}
	}


	

}



void Grid::equations() {

	//pêtla po czasie, w której rozwi¹zuje uk³ad równañ
	for (int tau = 0; tau < d.get_sim_time(); tau += d.get_sim_step()) {
		for (int i = 0; i < d.get_nN(); i++) {
			for (int j = 0; j < d.get_nN(); j++) {
				H_matrix_agr[i][j] = 0;
				C_matrix_agr[i][j] = 0;

			}
			//wpisanie temperatur z wêz³ów do wektora t0
			t0_vect[i] = nd[i].get_t();
			P_vect_agr[i] = 0;
		}
		//tworzenie zaagregowanych macierzy H,C oraz wektor P
		create_agregation();


		//wpisanie do macierzy H=H+C/dt
		for (int i = 0; i < d.get_nN(); i++) {
			for (int j = 0; j < d.get_nN(); j++) {
				H_matrix_agr[i][j] += (C_matrix_agr[i][j] / d.get_sim_step());
			}
		}

		//wpisanie do wektora P=C/dt-P
		for (int i = 0; i < d.get_nN(); i++) {
			double s = 0;
			for (int j = 0; j < d.get_nN(); j++) {
				s += (C_matrix_agr[i][j] / d.get_sim_step()) * t0_vect[j];
			}
			P_vect_agr[i] = s - P_vect_agr[i];
		}

		//stworzenie macierzy pomocniczej, która ma nN wierszy oraz nN+1 kolumn
		//wpisje do niej wartoœci z macierzy H, a w ostatniej kolumnie wartoœci z wektora P
		for (int i = 0; i < d.get_nN(); i++) {
			for (int j = 0; j < d.get_nN(); j++) {
				mat[i][j] = H_matrix_agr[i][j];
			}
		}
		for (int i = 0; i < d.get_nN(); i++) {
			mat[i][d.get_nN()] = P_vect_agr[i];
		}

		//wywo³anie funkcj rozwi¹zuj¹cej uk³ad równañ
		t1_vect = gaussianElimination(mat);

		cout << "Time: " << tau + d.get_sim_step();
		print_min_max_temp();




		//przypisanie rozwi¹zañ do poszczególnych wêz³ów
		for (int i = 0; i < d.get_nN(); i++) {
			nd[i].set_t(t1_vect[i]);
		}

	}
}
	



//funkcja do rowi¹zaywania uk³adu równañ metod¹ elimniacji Gaussa
//funkcje swap_row, forwardElim, backSub wszystkie s³u¿a do rozwi¹zywania uk³¹du równañ
double* Grid::gaussianElimination(double** mat)
{
	
	int singular_flag = forwardElim(mat);


	if (singular_flag != -1)
	{
		printf("Singular Matrix.\n");
		if (mat[singular_flag][d.get_nN()])
			printf("Inconsistent System.");
		else
			printf("May have infinitely many "
				"solutions.");

		return 0;
	}
	double* x = new double[d.get_nN()];
	
	x = backSub(mat);
	return x;
}


void Grid::swap_row(double** mat, int i, int j)
{


	for (int k = 0; k <= d.get_nN(); k++)
	{
		double temp = mat[i][k];
		mat[i][k] = mat[j][k];
		mat[j][k] = temp;
	}
}




int Grid::forwardElim(double** mat)
{
	for (int k = 0; k < d.get_nN(); k++)
	{
		
		int i_max = k;
		int v_max = abs(mat[i_max][k]);

		for (int i = k + 1; i < d.get_nN(); i++)
			if (fabs(mat[i][k]) > v_max)
				v_max = mat[i][k], i_max = i;

		if (!mat[k][i_max])
			return k; 

		
		if (i_max != k)
			swap_row(mat, k, i_max);


		for (int i = k + 1; i < d.get_nN(); i++)
		{
			double f = mat[i][k] / mat[k][k];

			for (int j = k + 1; j <= d.get_nN(); j++)
				mat[i][j] -= mat[k][j] * f;

			
			mat[i][k] = 0;
		}

		
	}
	
	return -1;
}


double* Grid::backSub(double** mat)
{
	
	double* x = new double[d.get_nN()];
	
	for (int i = d.get_nN() - 1; i >= 0; i--)
	{
		
		x[i] = mat[i][d.get_nN()];

		for (int j = i + 1; j < d.get_nN(); j++)
		{
		
			x[i] -= mat[i][j] * x[j];
		}

		x[i] = x[i] / mat[i][i];
	}
	return x;
	
}


//wypisanie elemntów oraz wêz³ów siatki
void Grid::print_grid() {
	
	
	cout << "Nodes" << endl;
	for (int i = 0; i < d.get_nN(); i++) {
		cout <<"Id: "<< i << " x:" << nd[i].get_x() << " y:" << nd[i].get_y() << " BC" << nd[i].get_BC() << endl;
	}

	cout << "\nElements\n";
	int a, b, c, e;
	for (int i = 0; i < d.get_nE(); i++) {
		
			a = elem[i].getID0();
			b = elem[i].getID1();
			c = elem[i].getID2();
			e = elem[i].getID3();

			cout << "\nElement "<<i;
			cout << endl;
			cout << " ID:" << a<< "- x:" << nd[a].get_x() << "y:" << nd[a].get_y() << "BC: "<<nd[a].get_BC()<<endl;
			cout << " ID: " << b<< "- x:" << nd[b].get_x() << "y:" << nd[b].get_y() << "BC: " << nd[b].get_BC() << endl;
			cout << " ID: " << c<< "- x:" << nd[c].get_x() << "y:" << nd[c].get_y() << "BC: " << nd[c].get_BC() << endl;
			cout << " ID: " <<e<< "- x:" << nd[e].get_x() << "y:" << nd[e].get_y() << "BC: " << nd[e].get_BC() << endl;


		


	}
}

//wypisywanie macierzy H, C oraz wektora P
void Grid::print_matrix() {
	cout << endl;
	cout << "matrix H :" << endl;
	for (int i = 0; i < d.get_nN(); i++) {
		for (int j = 0; j < d.get_nN(); j++) {
			cout << H_matrix_agr[i][j] << " ";
		}
		cout << endl;
	}
	
	cout << "matrix C :" << endl;
	for (int i = 0; i < d.get_nN(); i++) {
		for (int j = 0; j < d.get_nN(); j++) {
			cout << C_matrix_agr[i][j] << " ";
		}
		cout << endl;
	}

	
	cout << "Vector P:" << endl;
	for (int i = 0; i < d.get_nN(); i++) {
		cout << P_vect_agr[i]<<" ";
	}
	cout << endl;


}

//funkcja wyszukuj¹ca min i max temperatury oraz wypisuj¹ca je
void Grid::print_min_max_temp() {
	double min = t1_vect[0];

	for (int i = 1; i < d.get_nN(); i++) 
		if (min > t1_vect[i])
			min = t1_vect[i];
	double max = t1_vect[0];
	for (int i = 1; i < d.get_nN(); i++) 
		if (max < t1_vect[i])
			max = t1_vect[i];

	cout << "      Min Temp:  " << min << "   Max Temp:  " << max << endl;
}






Grid::~Grid() {
	for (int i = 0; i < d.get_nN(); i++) {
		delete[]H_matrix_agr[i];
		delete[]C_matrix_agr[i];
	}
	delete[]H_matrix_agr;
	delete[]C_matrix_agr;
	delete[]nd;
	delete[]elem;
	delete[]l;
}




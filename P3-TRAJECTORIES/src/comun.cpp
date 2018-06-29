/* 
 * File:   comun.cpp
 * Author: juane
 *
 * Created on 19 May 2018
 * 
 * USO: ./program <"name.dat"> <seed>
 */

/*INCLUDES*/

#include <fstream>
#include <algorithm>
#include <chrono> // std::chrono::system_clock
#include <dirent.h>
#include <cmath>

#include "comun.h"
#include "timer.h"
#include "random.h"


/////////////////////END INCLUDES////////////////////

//
//int evs_made = 0;

/*Read matrices from file "filep"
 * Suposse that format of file is:
 *	n -> size of square matrix
 *	m1
 *	m2
 */
bool readMatrices(vector<vector<int> > &f, vector<vector<int> > &d, string filep) {
	int size_p;
	int aux;
	string line;
	ifstream file(filep);

	if (file.is_open()) {
		f.clear();
		d.clear();
		file >> size_p;
		f.resize(size_p);
		d.resize(size_p);

		while (!file.eof()) {
			for (int i = 0; i < size_p; i++) {
				for (int j = 0; j < size_p; j++) {
					file>>aux;
					f[i].push_back(aux);
				}
			}

			for (int i = 0; i < size_p; i++) {
				for (int j = 0; j < size_p; j++) {
					file>>aux;
					d[i].push_back(aux);
				}
			}
		}

		file.close();
		return true;
	} else
		return false;
}

/////////////*END DEFINES, CONSTANTS AND GLOBAL FUNCTIONS*/////////////////////

/*Return cost of solution S*/
int cost(const vector<int> &S, const vector<vector<int> > &f, const vector<vector<int> > &d) {
	int cost = 0;

	//the size of f is equal that d
	for (int i = 0; i < f.size(); i++) {
		for (int j = 0; j < f.size(); j++) {
			if (j != i) {
				cost += f[i][j] * d[S[i]][S[j]];
			}
		}
	}

	return cost;
}


//Reset an fill the vector with sequence from 1 to nelem

void fillVector(vector<int> &v, int tam) {
	v.resize(tam);

	for (int i = 0; i < v.size(); ++i) {
		v[i] = i;
	}

	std::random_shuffle(v.begin(), v.end()); //mezclamos solucion inicial 

	//	for (int i = 0; i < v.size(); i++) {
	//		cout << v[i] << " ";
	//	}
	//	cout << endl ;
	//	cin.get();
}


///////////////////////UTILS BL/////////////////////////////////

/*Cost-Fact: cost of swap r and s in S*/
int costFact(const vector<int> &S, const vector<vector<int> > &f, const vector<vector<int> > &d, int r, int s) {
	int cost = 0;

	//the size of f is equal that d
	for (int k = 0; k < f.size(); k++) {
		if (k != r && k != s) {
			cost +=
					f[r][k] * (d[S[s]][S[k]] - d[S[r]][S[k]]) +
					f[s][k] * (d[S[r]][S[k]] - d[S[s]][S[k]]) +
					f[k][r] * (d[S[k]][S[s]] - d[S[k]][S[r]]) +
					f[k][s] * (d[S[k]][S[r]] - d[S[k]][S[s]]);
		}
	}

	return cost;
}

/*Operator neighbor: obtain a new neighbor*/
void Neigh(vector<int> &S, int i, int j) {
	int aux = S[i];
	S[i] = S[j];
	S[j] = aux;
}

/*Local Search*/
Solucion BL(const Solucion &ini, const vector<vector<int> > &f, const vector<vector<int> > &d) {
	Solucion SOL;
	int n_evals = 0;
	vector<int> best_perm(ini.first);

	int best_cost = cost(best_perm, f, d);

	vector<bool> mask(best_perm.size(), false);
	bool improve_flag = true;

	while (improve_flag && n_evals < MAX_EVALS) {
		for (int i = 0; i < ini.first.size() && n_evals < MAX_EVALS; ++i) {//Permutaciones posibles --> n*(n-1)/2
			if (!mask[i]) {
				improve_flag = false;

				for (int j = 0; j < ini.first.size() && n_evals < MAX_EVALS; ++j) {
					if (j != i) {
						//check move
						int cost_neigh = costFact(best_perm, f, d, i, j);
						n_evals++;

						if (cost_neigh < 0) { //if improves
							Neigh(best_perm, i, j);
							best_cost += cost_neigh;
							mask[i] = false;
							mask[j] = false;

							improve_flag = true;
						}

						//end check move
						//cout << n_evals << endl;
					}
				}
				if (!improve_flag)
					mask[i] = true;
			}
		}
	}

	//cout << n_evals << endl;
	SOL.first = best_perm;
	SOL.second = best_cost;

	return SOL;
}

////////////////////*END UTILS BL*//////////////////////////////

//Operators

/* Mut operator: swap gen with other random pos */
void mutOp(Solucion& S, int subl) {
	int rand_ini = Randint(0, subl - 1);
	int rand_fin = (rand_ini + subl) % (S.first.size() - 1);

	random_shuffle(S.first.begin() + rand_ini, S.first.begin() + rand_fin);
}

//GREEDY UTILS//

/*Return ^fi, i: 0..n-1*/
int fi(int i, const vector<vector<int> > &m) {
	int fi = 0;
	for (int j = 0; j < m.size(); j++)
		fi += m[i][j] + m[j][i];

	return fi;
}

/*Return ^dk, k: 0..n-1*/
int dk(int i, const vector<vector<int> > &m) {
	return fi(i, m);
}

/*Greedy*/
Solucion Greedy(const vector<vector<int> > &f, const vector<vector<int> > &d) {
	Solucion SOL;
	int size_p = f.size(), curr_fi, curr_dk;
	int max_fi = -1, min_dk = MAX_INT, max_fi_ind, min_dk_ind;
	vector<int> pf, pd; //potentials of d and f
	vector<bool> chpf, chpd; //cached potentials f and d already calculated

	vector<int> S; //vector of SOLUTIONS: ui->lk where s[i]=lk

	S.resize(size_p, 0);
	chpd.resize(size_p, false);
	chpf.resize(size_p, false);
	pf.resize(size_p, 0);
	pd.resize(size_p, 0);

	//GREEDY

	//calc potentials fi and dk and first solution
	for (int i = 0; i < size_p; i++) {
		curr_fi = fi(i, f);
		curr_dk = dk(i, d);

		//store current calcs in pf and pd
		pf[i] = curr_fi;
		pd[i] = curr_dk;

		if (curr_fi > max_fi) {
			max_fi_ind = i;
			max_fi = curr_fi;
		}

		if (curr_dk < min_dk) {
			min_dk_ind = i;
			min_dk = curr_dk;
		}
	}

	//store first solution
	S[max_fi_ind] = min_dk_ind;

	//caching fi and dk
	chpf[max_fi_ind] = true;
	chpd[min_dk_ind] = true;

	//from 1 to n-1(first solution already calc) calc remaining solutions and store
	for (int i = 0; i < size_p - 1; i++) {
		max_fi = -1, min_dk = MAX_INT;

		//the size of pf and pd is equal, calc greater fi and less dk
		for (int j = 0; j < pf.size(); j++) {
			if (!chpf[j]) {
				curr_fi = pf[j];
				if (curr_fi > max_fi) {
					max_fi_ind = j;
					max_fi = curr_fi;
				}
			}
			if (!chpd[j]) {
				curr_dk = pd[j];
				if (curr_dk < min_dk) {
					min_dk_ind = j;
					min_dk = curr_dk;
				}
			}
		}
		S[max_fi_ind] = min_dk_ind;
		chpf[max_fi_ind] = true;
		chpd[min_dk_ind] = true;
	}

	//END GREEDY

	SOL.first = S;
	SOL.second = cost(S, f, d);

	return SOL;
}

int cost_add(const vector<int> &S, const vector<vector<int> > &f, const vector<vector<int> > &d, int i, int k) {
	int cost = 0;

	for (int j = 0; j < S.size(); j++) {
		if (S[j] != -1)
			cost += f[i][S[j]] * d[k][S[j]];
	}

	//	if(cost==0)
	//		cost=0x0fffffff;

	return cost;
}

bool inS(const vector<int> &S, const assign &as) {
	if (S[as.u] != -1)
		return true;
	for (int i = 0; i < S.size(); i++) {
		if (S[i] == as.l)
			return true;
	}


	return false;
}
//END GREEDY UTILS

//////////////////////////*ALGORITHMS*////////////////////////////

//Simulated Annealing

Solucion ES(const Solucion &ini, vector<vector<int> > f, vector<vector<int> > d) {
	Solucion sol, best_sol;
	int n_evals = 0;

	int tam = f.size();
	int enfriamientos = 0;
	double crit_accept;

	double mu = 0.3, theta = 0.3;

	sol = ini;
	sol.second = cost(sol.first, f, d);

	double temp_ini = sol.second * mu / -log(theta);
	double temp_fin = 1E-3;

	int max_vecinos = 10 * tam;
	int max_exitos = 0.1 * max_vecinos;
	int M = MAX_EVALS / max_vecinos;
	double beta = (temp_ini - temp_fin) / (M * temp_ini * temp_fin);

	best_sol = sol;
	bool seguir_enfriando = true;
	double temp_act = temp_ini;

	//
	//	/*Una buena señal para el funcionamiento del algoritmo seria que beta fuera cercano a cero, 
	//	 * con estas condiciones, donde M es mayor a 4, la temperatura desciende rapidisimo, por lo que
	//	 * el algoritmo no funciona bien, deja de aceptar peores soluciones de manera muy rapida, por
	//	 * lo que decido utilizar el otro esquema de enfriamiento
	//	 */
	//

	int exitos = 0;
	int neighs_gen = 0;
	bool seguir_generando = true;

	while (seguir_enfriando) {

		seguir_generando = true;

		neighs_gen = 0;
		exitos = 0;

		while (seguir_generando) {
			neighs_gen++;

			int r, s;
			r = Randint(0, tam - 1);
			do {
				s = Randint(0, tam - 1);
			} while (r == s);

			long cost_fact = costFact(sol.first, f, d, r, s);
			n_evals++;

			crit_accept = exp(-cost_fact / temp_act);

			if (cost_fact < 0 || Randfloat(0, 1) <= crit_accept) {
				Neigh(sol.first, r, s);

				sol.second += cost_fact;

				if (sol.second < best_sol.second) {
					best_sol = sol;
				}

				exitos++;
			}

			if (neighs_gen >= max_vecinos || exitos >= max_exitos) {
				seguir_generando = false;

			}

		}

		if (n_evals >= MAX_EVALS || !exitos) { // || temp_act <= temp_fin
			seguir_enfriando = false;
			//cout << "SE PARA\n";
		} else {
			enfriamientos++;

			double alfa = Randfloat(0.9, 0.99);
			temp_act = temp_act * alfa;
			//temp_act = temp_act / (1 + beta * temp_act);
		}

	}

	//cout << "En " << enfriamientos << " enfriamientos. n_evals: " << n_evals << " .\n";
	return best_sol;
}

//GRASP

Solucion GRASP(const vector<vector<int> > &f, const vector<vector<int> > &d) {
	Solucion best_sol;

	long size_p, curr_fi, curr_dk;
	long max_fi = -1, min_fi = MAX_INT, min_dk = MAX_INT, max_dk = -1, max_fi_ind, min_dk_ind;

	vector<unsigned long> pf, pd; //potentials of d and f
	vector<int> lrcu, lrcl;

	double umbral_u, umbral_l, umbral_et2;
	double alfa = 0.3;

	Solucion partial_S; //vector of SOLUTIONS: Si->lk where S[i]=lk

	size_p = f.size();
	partial_S.first.resize(size_p, -1);
	pf.resize(size_p, 0);
	pd.resize(size_p, 0);

	//calc potentials fi and dk and first solution
	for (int i = 0; i < size_p; i++) {
		curr_fi = fi(i, f);
		curr_dk = dk(i, d);

		//store current calcs in pf and pd
		pf[i] = curr_fi;
		pd[i] = curr_dk;


		if (curr_fi > max_fi) {
			max_fi_ind = i;
			max_fi = curr_fi;
		}
		if (curr_fi < min_fi) {
			min_fi = curr_fi;
		}

		if (curr_dk < min_dk) {
			min_dk_ind = i;
			min_dk = curr_dk;
		}
		if (curr_dk > max_dk)
			max_dk = curr_dk;
	}


	umbral_u = max_fi - alfa * (max_fi - min_fi);
	umbral_l = min_dk + alfa * (max_dk - min_dk);

	for (int i = 0; i < size_p; i++) {
		if (pf[i] >= umbral_u)
			lrcu.push_back(i);
		if (pd[i] <= umbral_l)
			lrcl.push_back(i);
	}

	int size_lrcu = lrcu.size(), size_lrcl = lrcl.size();

	//ETAPA1
	int cu1, cu2, cl1, cl2;

	for (int w = 0; w < MAX_ITERS; w++) {
		//candidatos lrcu y lrcl: puede ser que, debido a los umbrales, no dispongamos de dos
		//candidatos aleatorios a coger.

		if (size_lrcu > 1 && size_lrcl > 1) {
			cu1 = Randint(0, size_lrcu - 1);
			do {
				cu2 = Randint(0, size_lrcu - 1);
			} while (cu2 == cu1);

			cl1 = Randint(0, size_lrcl - 1);
			do {
				cl2 = Randint(0, size_lrcl - 1);
			} while (cl2 == cl1);

			partial_S.first[lrcu[cu1]] = lrcl[cl1];
			partial_S.first[lrcu[cu2]] = lrcl[cl2];
		} else {
			if (size_lrcu == 1 && size_lrcl > 1) {
				cu1 = 0;

				cl1 = Randint(0, size_lrcl - 1);
				do {
					cl2 = Randint(0, size_lrcl - 1);
				} while (cl2 == cl1);

				do {
					cu2 = Randint(0, size_p - 1);
				} while (cu2 == lrcu[cu1]);

				partial_S.first[lrcu[cu1]] = lrcl[cl1];
				partial_S.first[cu2] = lrcl[cl2];
			} else if (size_lrcu > 1 && size_lrcl == 1) {
				cu1 = Randint(0, size_lrcu - 1);
				do {
					cu2 = Randint(0, size_lrcu - 1);
				} while (cu2 == cu1);

				cl1 = 0;
				do {
					cl2 = Randint(0, size_p - 1);
				} while (cl2 == lrcl[cl1]);

				partial_S.first[lrcu[cu1]] = lrcl[cl1];
				partial_S.first[lrcu[cu2]] = cl2;
			} else {
				cout << "ERROR!\n\n";
				break;
			}
		}


		//FIN ETAPA1

		//ETAPA2

		for (int u = 0; u < size_p - 2; u++) {
			vector<assign> lc2;

			for (int i = 0; i < size_p; i++) {
				for (int k = 0; k < size_p; k++) {
					assign aux;
					aux.u = i;
					aux.l = k;

					if (!inS(partial_S.first, aux)) {
						aux.cost_add = cost_add(partial_S.first, f, d, i, k);

						lc2.push_back(aux);
					}
				}
			}

			sort(lc2.begin(), lc2.end());

			int min_lc2 = lc2[0].cost_add, max_lc2 = lc2[lc2.size() - 1].cost_add;

			umbral_et2 = min_lc2 + alfa * (max_lc2 - min_lc2);

			vector<assign> lrc2;

			for (int i = 0; i < lc2.size(); i++) {
				if (lc2[i].cost_add <= umbral_et2)
					lrc2.push_back(lc2[i]);
				else
					break;
			}

			int rand_lrc2 = Randint(0, lrc2.size() - 1);

			partial_S.first[lrc2[rand_lrc2].u] = lrc2[rand_lrc2].l;
		}


		Solucion act_sol = BL(partial_S, f, d);

		if (act_sol.second < best_sol.second)
			best_sol = act_sol;
		//FIN ETAPA2

		partial_S.first.clear();
		partial_S.first.resize(size_p, -1);
	}

	return best_sol;
}

//ILS

Solucion ILS(const vector<vector<int> > &f, const vector<vector<int> > &d) {
	Solucion sol_ini, best_sol;

	int size_p = f.size();
	int subl = size_p / 4;

	fillVector(sol_ini.first, size_p);

	sol_ini = BL(sol_ini, f, d);

	for (int i = 0; i < MAX_ITERS - 1; i++) {
		if (sol_ini.second < best_sol.second)
			best_sol = sol_ini;

		mutOp(best_sol, subl);
		sol_ini = BL(best_sol, f, d);
	}

	return best_sol;
}

//BMB

Solucion BMB(const vector<vector<int> > &f, const vector<vector<int> > &d) {
	int size_p = f.size();
	Solucion best_sol;

	for (int i = 0; i < MAX_ITERS; i++) {
		Solucion s_act;
		fillVector(s_act.first, size_p);

		s_act = BL(s_act, f, d);

		if (s_act.second < best_sol.second)
			best_sol = s_act;
	}

	return best_sol;
}

//ILS-ES

Solucion ILS_ES(const vector<vector<int> > &f, const vector<vector<int> > &d) {
	Solucion sol_ini, best_sol;

	int size_p = f.size();
	int subl = size_p / 4;

	fillVector(sol_ini.first, size_p);
	sol_ini = ES(sol_ini, f, d);

	for (int i = 0; i < MAX_ITERS - 1; i++) {
		if (sol_ini.second < best_sol.second)
			best_sol = sol_ini;

		mutOp(best_sol, subl);

		sol_ini = ES(best_sol, f, d);
	}

	return best_sol;
}


///////////////////////*END ALGORITHMS*///////////////////////////



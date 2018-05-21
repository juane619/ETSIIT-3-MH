/* 
 * File:   main.cpp
 * Author: juane
 *
 * Created on 07 March 2018, 23:44
 * DOING AGG
 * USO: ./program <"name.dat"> <seed>
 */

/*INCLUDES*/

#include <fstream>
#include <algorithm>
#include <chrono> // std::chrono::system_clock
#include <dirent.h>

#include "comun.h"
#include "timer.h"
#include "random.h"


/////////////////////END INCLUDES////////////////////
int evs_made = 0;

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

/*GENETICS ALGORITHMS*/

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
genoma BL(const vector<int> &ini, const vector<vector<int> > &f, const vector<vector<int> > &d) {
	genoma SOL;
	int n_evals = 0;
	vector<int> best_perm(ini);

	int best_cost = cost(best_perm, f, d);
	evs_made++;

	vector<bool> mask(best_perm.size(), false);
	bool improve_flag = true;

	while (improve_flag && n_evals < MAX_EVAL_NEIGH) {
		for (int i = 0; i < ini.size(); ++i) {//Permutaciones posibles --> n*(n-1)/2
			if (!mask[i]) {
				improve_flag = false;

				for (int j = 0; j < ini.size(); ++j) {
					if (j != i) {
						//check move
						int cost_neigh = costFact(best_perm, f, d, i, j);
						//evs_made++;

						if (cost_neigh < 0) { //if improves
							Neigh(best_perm, i, j);
							best_cost += cost_neigh;
							mask[i] = false;
							mask[j] = false;

							improve_flag = true;
						}

						//end check move
						n_evals++;
						//cout << n_evals << endl;
					}
				}
				if (!improve_flag)
					mask[i] = true;
			}
		}
	}
	SOL.first = best_perm;
	SOL.second = best_cost;

	return SOL;
}

////////////////////*END UTILS BL*//////////////////////////////


//other functions
//get best crom

void getBest(const vector<genoma>& pop, genoma &best_gen, int &best) {
	int bestl = 0xffffffff;

	for (int i = 0; i < pop.size(); i++) {
		if (pop[i].second < bestl) {
			bestl = pop[i].second;
			best_gen = pop[i];
			best = i;
		}
	}
}

//selection by tournament

int tournament(int k, const vector<genoma>& pop_curr, vector<int>& parents) {
	int index_best = -1, cost_best = 0xffffffff;
	int curr_par, old_par = -1;
	//cout << "Lanzando torneo: " << K_TOUR << endl;

	for (int i = 0; i < k; i++) {
		do {
			curr_par = Randint(0, pop_curr.size() - 1);
		} while (curr_par == old_par);
		old_par = curr_par;

		//cout << pop_curr[curr_par].second << endl;

		if (pop_curr[curr_par].second < cost_best) {
			cost_best = pop_curr[curr_par].second;
			index_best = curr_par;
		}
	}

	//	cout << "ganador: " << index_best << " con " << cost_best << " puntos.\n";
	//	cin.get();
	return index_best;
}

/* preserve elitism: replace the best by the worst in new population 
 * precondition: best genoma before calculate
 */
void preserveElitism(const genoma& best, vector<genoma>& new_pop) {
	auto it = max_element(new_pop.begin(), new_pop.end());

	it->first = best.first;
	it->second = best.second;
}

//end other functions

//evaluate pop

void evalPop(vector<genoma>& pop, const vector<vector<int> > &f, const vector<vector<int> > &d) {
	for (int i = 0; i < pop.size(); i++) {
		if (pop[i].second == 0) {
			//cout << "evaluando..\n";
			pop[i].second = cost(pop[i].first, f, d);
		}

		evs_made++;
	}
}

//Start population

void generateStartPop(int size_gen, vector<genoma>& pop, int &best, const vector<vector<int> > &f, const vector<vector<int> > &d) {
	genoma gen;
	int best_cost = 0xffffffff, curr_cost;
	best = 0;

	for (int i = 0; i < SIZE_POP; i++) {
		fillVector(gen.first, size_gen);
		curr_cost = cost(gen.first, f, d);

		if (curr_cost < best_cost) {
			best_cost = curr_cost;
			best = i;
		}
		gen.second = curr_cost;
		pop.push_back(gen);

		evs_made++;
	}
}


//Selectg

void selectg(const vector<genoma>& pop_curr, vector<int>& parents) {
	for (int i = 0; i < pop_curr.size(); i++) {
		parents.push_back(tournament(K_TOUR, pop_curr, parents));
	}
}

//Selects

void selects(const vector<genoma>& pop_curr, vector<int>& parents) {
	for (int i = 0; i < 2; i++) {
		parents.push_back(tournament(K_TOUR, pop_curr, parents));
	}
}

//Operators

/* Mut operator: swap gen with other random pos */
void mutOp(genoma& S, const vector<vector<int> > &f, const vector<vector<int> > &d) {
	int sample_other, sample_gen;

	sample_gen = Randint(0, S.first.size() - 1);

	do {
		sample_other = Randint(0, S.first.size() - 1);
	} while (sample_other == sample_gen);

	//mutation
	Neigh(S.first, sample_gen, sample_other);

	S.second -= costFact(S.first, f, d, sample_gen, sample_other);
}

/* Cross operator: based-position
 * return true if best is replaced, else false
 */
bool crossOpPos(const genoma &par1, const genoma &par2, vector<genoma>& offsprings) {
	bool crossed = false;

	genoma offspring1, offspring2;
	vector<int> rest;

	int size_gen = par1.first.size();
	offspring1.first.resize(size_gen, -1);
	offspring2.first.resize(size_gen, -1);

	//true cross
	for (int k = 0; k < size_gen; k++) {
		if (par1.first[k] == par2.first[k]) {
			offspring1.first[k] = par1.first[k];
			offspring2.first[k] = par1.first[k];
			crossed = true;
		} else
			rest.push_back(par1.first[k]);
	}

	if (!rest.empty() && crossed) {
		int sz_rest = rest.size();
		random_shuffle(rest.begin(), rest.end());

		int k = 0;
		for (int l = 0; l < size_gen && k < sz_rest; l++) {
			if (offspring1.first[l] == -1)
				offspring1.first[l] = rest[k++];
		}

		offsprings.push_back(offspring1);

		//metemos uno de los dos padres al azar
		//		if (Randint(1, 2) == 1)
		//			offsprings.push_back(par1);
		//		else
		//			offsprings.push_back(par2);

		//metemos al mejor padre
		if (par1.second < par2.second)
			offsprings.push_back(par1);
		else
			offsprings.push_back(par2);

		//metemos otro hijo permutando
		//		random_shuffle(rest.begin(), rest.end());
		//		k = 0;
		//		for (int l = 0; l < size_gen && k < sz_rest; l++) {
		//			if (offspring2.first[l] == -1)
		//				offspring2.first[l] = rest[k++];
		//		}
		//
		//		offsprings.push_back(offspring2);

		return true;
	} else { //if no cross, put two croms
		offsprings.push_back(par1);
		offsprings.push_back(par2);
		return false;
	}
}

/* Cross operator:PMX
 * return true if best is replaced, else false
 */
bool crossOpPMX(const genoma &par1, const genoma &par2, vector<genoma>& offsprings) {
	if (par1.first != par2.first) {
		int cutpoint1, cutpoint2;
		genoma offspring1, offspring2;
		vector<int> sub1, sub2;

		int size_gen = par1.first.size();

		offspring1.first.resize(size_gen, -1);
		offspring2.first.resize(size_gen, -1);

		cutpoint1 = Randint(0, size_gen - 1);
		do {
			cutpoint2 = Randint(0, size_gen - 1);
		} while (cutpoint1 == cutpoint2);

		if (cutpoint1 > cutpoint2) {
			int aux = cutpoint1;
			cutpoint1 = cutpoint2;
			cutpoint2 = aux;
		}

		//true cross
		//create substrings and start fill offsprings
		for (int k = cutpoint1; k <= cutpoint2; k++) { //sure not outside bounds of vector
			sub1.push_back(par1.first[k]);
			offspring1.first[k] = par2.first[k];
			sub2.push_back(par2.first[k]);
			offspring2.first[k] = par1.first[k];
		}
		//end create substrings and start fill offsprings

		//cross offs1 and offs2
		for (int k = 0; k < size_gen; k++) { //sure not outside bounds of vector
			if (k < cutpoint1 || k > cutpoint2) {
				vector<int>::iterator it;

				if ((it = find(sub2.begin(), sub2.end(), par1.first[k])) == sub2.end())
					offspring1.first[k] = par1.first[k];
				else {
					int aux_index, new_k;

					do {
						aux_index = it - sub2.begin();
						new_k = sub1[aux_index];
					} while ((it = find(sub2.begin(), sub2.end(), new_k)) != sub2.end());

					offspring1.first[k] = new_k;
				}
				if ((it = find(sub1.begin(), sub1.end(), par2.first[k])) == sub1.end())
					offspring2.first[k] = par2.first[k];
				else {
					int aux_index, new_k;

					do {
						aux_index = it - sub1.begin();
						new_k = sub2[aux_index];
					} while ((it = find(sub1.begin(), sub1.end(), new_k)) != sub1.end());

					offspring2.first[k] = new_k;
				}
			}
		}

		offsprings.push_back(offspring1);
		offsprings.push_back(offspring2);

		return true;
	} else {
		offsprings.push_back(par1);
		offsprings.push_back(par2);
		return false;
	}
}

/*run genetic generational algorithm*/
genoma AGG(const vector<vector<int> > &f, const vector<vector<int> > &d, int size_gen, int TYPE_AGG) {
	bool parada = false, lost_cross = true, lost_mut = false;
	int cont = 0, old_best_cost = 0xffffffff;
	vector<genoma> curr_pop;
	int best = -1;

	int n_cross = PROB_CROSS_AGG * SIZE_POP / 2;
	int n_mut = PROB_MUT * SIZE_POP*size_gen;

	//	cout << n_cross << " " << n_mut << " " << evs_made << endl;
	//	cin.get();

	//initial population
	generateStartPop(size_gen, curr_pop, best, f, d);
	genoma best_gen = curr_pop[best];

	int n_generations = 0;

	while (!parada) {
		//select parents to cross
		vector<int> parents;
		selectg(curr_pop, parents);

		//cross and replace
		int i = 0, j = 0;
		vector<genoma> new_pop;
		lost_cross = true;

		while (i < parents.size()) {
			genoma par1, par2;
			if (j < n_cross) {
				par1 = curr_pop[parents[i]], par2 = curr_pop[parents[i + 1]];

				//cross
				switch (TYPE_AGG) {
					case 1:
						crossOpPos(par1, par2, new_pop);
						break;
					case 2:
						crossOpPMX(par1, par2, new_pop);
						break;
					default:
						cout << "TIPO INCORRECTO AGG.\n";
				}

				j++;
				i += 2;
			} else {
				if (parents[i] == best) {
					lost_cross = false;
				}
				new_pop.push_back(curr_pop[parents[i++]]);
			}
		}
		//end cross

		//eval only new additions
		evalPop(new_pop, f, d);

		//mutation
		int sample_crom;

		for (int i = 0; i < n_mut; i++) {
			//cout << "mutacion AGG: " << i << ": \n";
			sample_crom = Randint(0, new_pop.size() - 1);
			genoma S = new_pop[sample_crom];

			if (S == best_gen)
				lost_mut = true;

			mutOp(S, f, d);

			new_pop[sample_crom] = S;
		}
		//end mutation

		//if the best lost, replace by the worst of this population
		if (lost_mut || lost_cross) {
			preserveElitism(best_gen, new_pop);
			lost_mut = false;
		}

		getBest(new_pop, best_gen, best);

		//stop condition
		if (evs_made > MAX_EVALS)
			parada = true;
		else {
			curr_pop = new_pop;

			//			int dif = old_best_cost - best_gen.second;
			//
			//			if (dif < 0 && n_generations > 0) {
			//				cout << "FALLO en iteracion " << n_generations << endl << endl;
			//				for (int i = 0; i < curr_pop.size(); i++) {
			//					cout << curr_pop[i].second << endl;
			//				}
			//				cout << endl;
			//
			//				cin.get();
			//			}
			//
			//			if (dif == 0)
			//				cont++;
			//			else {
			//				cont = 0;
			//			}

			//			if (cont > GEN_NO_IMPROVE_AGG) { //GENERACIONES SIN MEJORAR
			//				parada = true;
			//			} else {
			//				old_best_cost = best_gen.second;
			//			}
		}

		//		cout << n_generations << " " << evs_made << " " << best_gen.second << endl;
		//		cin.get();
		n_generations++;
	}

	cout << "Se han hecho " << n_generations << " generaciones y " << evs_made << " evaluaciones.\n";
	//	cin.get();
	evs_made = 0;
	return best_gen;
}

/*stationary genetic*/
genoma AGS(const vector<vector<int> > &f, const vector<vector<int> > &d, int size_gen, int TYPE_AGG) {
	bool parada = false;
	int cont = 0, old_best_cost = 0xffffffff;
	vector<genoma> curr_pop;
	int best = 0;
	genoma best_gen;

	double prob_mut = PROB_MUT * 2 * size_gen;

	//initial population
	generateStartPop(size_gen, curr_pop, best, f, d);

	int n_generations = 0;
	while (!parada) {

		//select parents to cross
		vector<int> parents;
		selects(curr_pop, parents);

		//cross and replace
		genoma par1, par2;
		par1 = curr_pop[parents[0]], par2 = curr_pop[parents[1]];

		//cross
		vector<genoma> new_pop;

		//cross
		switch (TYPE_AGG) {
			case 1:
				crossOpPos(par1, par2, new_pop);
				break;
			case 2:
				crossOpPMX(par1, par2, new_pop);
				break;
			default:
				cout << "TIPO INCORRECTO AGG.\n";
		}

		evalPop(new_pop, f, d);

		//mutation
		int sample_crom;

		if (Rand() <= prob_mut) {
			sample_crom = Randint(0, new_pop.size() - 1);

			genoma S = new_pop[sample_crom];

			mutOp(S, f, d);
			new_pop[sample_crom] = S;
		}
		//end mutation

		sort(curr_pop.begin(), curr_pop.end());
		sort(new_pop.begin(), new_pop.end());

		//replace if improve
		if (new_pop[0].second < curr_pop[SIZE_POP - 2].second)
			curr_pop[SIZE_POP - 2] = new_pop[0];

		if (new_pop[1].second < curr_pop[SIZE_POP - 1].second)
			curr_pop[SIZE_POP - 1] = new_pop[1];

		//get best gen of population and index in population
		getBest(curr_pop, best_gen, best);

		//eval to follow or stop
		if (evs_made > MAX_EVALS)
			parada = true;
		else {
			//cout << endl << best_gen.second << endl;
			//			int dif = old_best_cost - best_gen.second;
			//
			//			if (dif < 0 && n_generations > 0) {
			//				cout << "FALLO en iteracion " << n_generations << endl << endl;
			//				//				for (int i = 0; i < curr_pop.size(); i++) {
			//				//					cout << curr_pop[i].second << endl;
			//				//				}
			//				//				cout << endl;
			//
			//				cin.get();
			//			}
			//
			//			if (dif == 0) {
			//				cont++;
			//			} else {
			//				cont = 0;
			//			}
			//
			//			if (cont > GEN_NO_IMPROVE_AGS) { //GENERACIONES SIN MEJORAR
			//				parada = true;
			//			} else
			//				old_best_cost = best_gen.second;

		}

		//cout << n_generations << " " << evs_made << " " << best_gen.second << endl;
		n_generations++;
	}
	cout << "Se han hecho " << n_generations << " generaciones y " << evs_made << " evaluaciones.\n";

	evs_made = 0;
	return best_gen;
}

/*run memetic algorithm*/
genoma AM(const vector<vector<int> > &f, const vector<vector<int> > &d, int size_gen, int TYPE_AGG, int TYPE_AM) {
	bool parada = false, lost_cross = true, lost_mut = false, mutando = false;
	int cont = 0, old_best_cost = 0xffffffff;
	vector<genoma> curr_pop;
	int best = -1;
	genoma best_gen;

	int n_cross = PROB_CROSS_AGG * SIZE_POP / 2;
	int n_mut = PROB_MUT * SIZE_POP*size_gen;

	//initial population
	generateStartPop(size_gen, curr_pop, best, f, d);
	best_gen = curr_pop[best];

	int n_generations = 0;

	while (!parada) {
		//select parents to cross
		vector<int> parents;
		selectg(curr_pop, parents);

		//cross and replace
		int i = 0, j = 0;
		vector<genoma> new_pop;
		lost_cross = true;

		while (i < parents.size()) {
			genoma par1, par2;
			if (j < n_cross) {
				par1 = curr_pop[parents[i]], par2 = curr_pop[parents[i + 1]];

				//cross
				switch (TYPE_AGG) {
					case 1:
						crossOpPos(par1, par2, new_pop);
						break;
					case 2:
						crossOpPMX(par1, par2, new_pop);
						break;
					default:
						cout << "TIPO INCORRECTO AGG.\n";
				}

				j++;
				i += 2;
			} else {
				if (parents[i] == best) {
					lost_cross = false;
				}
				new_pop.push_back(curr_pop[parents[i++]]);
			}
		}

		evalPop(new_pop, f, d);

		//end cross

		//mutation
		int sample_crom;

		for (int i = 0; i < n_mut; i++) {
			//cout << "mutacion AGG: " << i << ": \n";
			sample_crom = Randint(0, new_pop.size() - 1);
			genoma S = new_pop[sample_crom];

			if (S == best_gen)
				lost_mut = true;

			mutOp(S, f, d);

			new_pop[sample_crom] = S;
		}
		//end mutation

		//if the best lost, replace by the worst of this population
		if (lost_mut || lost_cross) {
			preserveElitism(best_gen, new_pop);
			lost_mut = false;
		}

		if ((n_generations + 1) % GENERATIONS_TO_MEM == 0) {
			//MEMETIC ALGORITHM
			vector<int> pop_mem;
			int size_pop_mem;

			switch (TYPE_AM) {
				case 1: //MEMETIC ALGORITHM (1)
					for (int k = 0; k < new_pop.size(); k++) {
						genoma new_bl = BL(new_pop[k].first, f, d);
						
						if (new_bl.second < new_pop[k].second)
							new_pop[k] = new_bl;
					}
					break;
				case 2:
					//MEMETIC ALGORITHM (2)
					for (int i = 0; i < new_pop.size(); i++) {
						if (Randfloat(0, 1) <= 0.1)
							pop_mem.push_back(i);
					}

					for (int k = 0; k < pop_mem.size(); k++) {
						genoma new_bl = BL(new_pop[pop_mem[k]].first, f, d);

						if (new_bl.second < new_pop[pop_mem[k]].second)
							new_pop[k] = new_bl;
					}
					break;
				case 3:
					//MEMETIC ALGORITHM (3)
					sort(new_pop.begin(), new_pop.end());

					size_pop_mem = PROB_POP_MEM * new_pop.size();

					for (int i = 0; i < size_pop_mem; i++) {
						pop_mem.push_back(i);
					}

					for (int k = 0; k < pop_mem.size(); k++) {
						genoma new_bl = BL(new_pop[pop_mem[k]].first, f, d);

						if (new_bl.second < new_pop[pop_mem[k]].second)
							new_pop[k] = new_bl;
					}

					break;
				default:
					cout << "TIPO MAL ELEGIDO\n";
			}

			//END MEMETIC ALGORITHM
		}

		//get best gen of population and index in population
		getBest(new_pop, best_gen, best);

		if (evs_made > MAX_EVALS)
			parada = true;
		else {
			curr_pop = new_pop;

			//			cout << "\nFINAL: el mejor es: " << best << " " << best_gen.second << endl;
			//			cin.get();

			//			int dif = old_best_cost - best_gen.second;
			//
			//			if (dif < 0 && n_generations > 0) {
			//				cout << "FALLO en iteracion " << n_generations << endl << endl;
			//			}
			//
			//			if (dif == 0)
			//				cont++;
			//			else {
			//				cont = 0;
			//			}
			//
			//			if (cont > 100) { //GENERACIONES SIN MEJORAR
			//				parada = true;
			//			} else
			//				old_best_cost = best_gen.second;
		}

		//cout << n_generations << " " << evs_made << " " << best_gen.second << endl;
		n_generations++;
	}
	cout << "Se han hecho " << n_generations << " generaciones y " << evs_made << " evaluaciones.\n";
	//cin.get();

	evs_made = 0;
	return best_gen;
}

///////////////////////*END GENETIC ALGORITHMS*///////////////////////////



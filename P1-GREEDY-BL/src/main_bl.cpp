/* 
 * File:   main.cpp
 * Author: juane
 *
 * Created on 07 March 2018, 23:44
 * GREEDY and BL ALREADY DONE!
 * USO: ./program <"name.dat"> <seed>
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <algorithm>
#include <random>
#include <chrono> // std::chrono::system_clock
#include "timer.h"

#define  pair_sol pair<vector<int>, unsigned long> 

using namespace std;

const int EVALS = 50000;
const int N_EXECS = 200;

/*STRUCT Measure to print the measurement of either execution
 *
 */
struct Measure {
    unsigned long best, worst;
    double mean_time;
    long long mean_cost;

    Measure() {
        worst = 0;
        best = 0xffffffff;
        mean_cost = 0;
        mean_time = 0.0;
    }

    void print() {
        cout << "Mejor caso: " << best << endl <<
                "Peor caso: " << worst << endl <<
                "Media coste: " << mean_cost << endl <<
                "Media tiempo: " << mean_time << endl;
    }
};

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
        return true;
    } else
        return false;
}


//UTILS BL

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

/*Local Search*/
pair_sol BL(const vector<int> &ini, const vector<vector<int> > &f, const vector<vector<int> > &d) {
    pair_sol SOL;
    int n_evals = 0;
    vector<int> best_perm(ini);

    int best_cost = cost(best_perm, f, d);

    vector<bool> mask(best_perm.size(), false);
    bool improve_flag = true;

    while (improve_flag && n_evals < EVALS) {
        for (int i = 0; i < ini.size(); ++i) {//Permutaciones posibles --> n*(n-1)/2
            if (!mask[i]) {
                improve_flag = false;

                for (int j = 0; j < ini.size(); ++j) {
                    if (j != i) {
                        //check move
                        int cost_neigh = costFact(best_perm, f, d, i, j);
                        if (cost_neigh < 0) { //if improves
                            Neigh(best_perm, i, j);
                            best_cost += cost_neigh;
                            mask[i] = false;
                            mask[j] = false;

                            improve_flag = true;
                        }
                        //end check move
                        n_evals++;
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

//Reset an fill the vector with sequence from 1 to nelem

void fillVector(vector<int> &v) {
    for (int i = 0; i < v.size(); ++i) {
        v[i] = i;
    }

    random_shuffle(v.begin(), v.end()); //mezclamos solucion inicial 
}

/*END UTILS BL*/

int main(int argc, char** argv) {
    if (argc < 3 || argc > 3) {
        cerr << "Uso: ./bl <file.data> <seed>" << endl;
        return -1;
    }
	
	int tam = 0;
    vector<vector<int> > f;
    vector<vector<int> > d;
    pair_sol sol_BL;
    vector<int> ini;
    double time_elapsed = -1;

    Measure med_bl = Measure();

    if (readMatrices(f, d, argv[1])) {
        //output file to boxplot
        string name(argv[1]);
        name = name.substr(11);
        string nameout("meds/");
        nameout += name;

        ofstream filout(nameout);

        tam = f.size();
        ini.resize(tam);
        unsigned seed = stoi(argv[2]);
        srand(seed);
        //EXEC BL

        // obtain a time-based seed:
        //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();	//Seed establecida a traves de la hora actual

        for (int i = 0; i < N_EXECS; i++) {
            fillVector(ini); //reset and fill the initial solution

            start_timers();
            sol_BL = BL(ini, f, d);
            time_elapsed = elapsed_time();

            if (sol_BL.second > med_bl.worst)
                med_bl.worst = sol_BL.second;
            if (sol_BL.second < med_bl.best)
                med_bl.best = sol_BL.second;
            med_bl.mean_time += time_elapsed;
            med_bl.mean_cost += sol_BL.second;

            //cout << sol_BL.second << endl;
            filout << sol_BL.second << endl;
        }
        filout.close();

        med_bl.mean_time = med_bl.mean_time / N_EXECS;
        med_bl.mean_cost = med_bl.mean_cost / N_EXECS;
        cout << endl;
        //PRINT S AND C(S) BL
        //
        //        for (int i = 0; i < sol_BL.first.size(); i++) {
        //            cout << sol_BL.first[i] << " ";
        //        }
        med_bl.print();
        //	cout << endl << sol_BL.second << endl;
        //END PRINT S AND C(S) BL
        //END EXEC BL

    } else
        cerr << "Can´t open the file\n";


    return 0;
}


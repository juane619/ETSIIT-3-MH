/* 
 * File:   main.cpp
 * Author: juane
 *
 * Created on 07 March 2018, 23:44
 * GREEDY and BL ALREADY DONE!
 * USO: ./program <"name.dat">
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <climits>
#include <algorithm>
#include <chrono> // std::chrono::system_clock
#include "timer.h"

#define  pair_sol pair<vector<int>, unsigned long> 

using namespace std;

const int EVALS = 50000;

/*STRUCT Measure to print the measurement of either execution
 *
 */
struct Measure {
    unsigned long best, worst;
    double mean;

    Measure() {
	worst = 0;
	best = 0xffffffff;
	mean = 0.0;
    }

    void print() {
	cout << "Mejor caso: " << best << endl <<
		"Peor caso: " << worst << endl <<
		"Media: " << mean << endl;
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

//UTILS GREEDY

/*Return cost of solution S*/
int cost(const vector<int> &S, const vector<vector<int> > &f, const vector<vector<int> > &d) {
    int cost = 0;

    //the size of f is equal that d
    for (int i = 0; i < f.size(); i++) {
	for (int j = 0; j < f.size(); j++) {
	    if (j != i) {
		cost += f[i][j] * d[S[i] - 1][S[j] - 1];
	    }
	}
    }
    return cost;
}

/*Return ^fi, i: 1..n*/
int fi(int i, const vector<vector<int> > &m) {
    int fi = 0;
    for (int j = 0; j < m.size(); j++) {
	fi += m[i - 1][j];
    }
    return fi;
}

/*Return ^dk, k: 1..n*/
int dk(int i, const vector<vector<int> > &m) {
    return fi(i, m);
}

/*Greedy-QAP*/
pair_sol Greedy(const vector<vector<int> > &f, const vector<vector<int> > &d) {
    pair_sol SOL;
    int size_p, curr_fi, curr_dk;
    int max_fi = -1, min_dk = INT_MAX, max_fi_ind, min_dk_ind;
    vector<int> pf, pd; //potentials of d and f
    vector<bool> chpf, chpd; //cached potentials f and d already calculated

    vector<int> S; //vector of SOLUTIONS: ui->lk where s[i]=lk

    size_p = f.size();
    S.resize(size_p, 0);
    chpd.resize(size_p, false);
    chpf.resize(size_p, false);
    pf.resize(size_p, 0);
    pd.resize(size_p, 0);

    //GREEDY

    //calc potentials fi and dk and first solution
    for (int i = 1; i <= size_p; i++) {
	curr_fi = fi(i, f);
	curr_dk = dk(i, d);

	//store current calcs in pf and pd
	pf[i - 1] = curr_fi;
	pd[i - 1] = curr_dk;

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
    S[max_fi_ind - 1] = min_dk_ind;

    //caching fi and dk
    chpf[max_fi_ind - 1] = true;
    chpd[min_dk_ind - 1] = true;

    //from 1 to n-1(first solution already calc) calc remaining solutions and store
    for (int i = 0; i < size_p - 1; i++) {
	max_fi = -1, min_dk = INT_MAX;

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
	S[max_fi_ind] = min_dk_ind+1;
	chpf[max_fi_ind] = true;
	chpd[min_dk_ind] = true;
    }

    //END GREEDY

    SOL.first = S;
    SOL.second = cost(S, f, d);

    return SOL;
}

//END UTILS GREEDY


int main(int argc, char** argv) {
    int tam = 0;
    vector<vector<int> > f;
    vector<vector<int> > d;
    pair_sol sol_greedy;
    vector<int> ini;
    double time_elapsed = -1;

    Measure med_gd = Measure();

    if (readMatrices(f, d, argv[1])) {
	tam = f.size();
	
	//EXEC GREEDY
	start_timers();
	sol_greedy = Greedy(f, d);
	time_elapsed = elapsed_time();
	med_gd.best = sol_greedy.second;
	med_gd.mean = time_elapsed;
	//PRINT S AND C(S) GREEDY
	med_gd.print();
	//		for (int i = 0; i < sol_greedy.first.size(); i++) {
	//			cout << sol_greedy.first[i] << " ";
	//		}
	//		cout << endl << sol_greedy.second << endl;
	//END PRINT S AND C(S) GREEDY

	//END EXEC GREEDY

    } else {
	cerr << "Can´t open the file\n";
    }

    return 0;
}


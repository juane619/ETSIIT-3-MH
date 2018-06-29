/* 
 * File:   main.cpp
 * Author: juane
 *
 * Created on 07 March 2018, 23:44
 * DOING AGG
 * USO: ./program <"name.dat"> <seed>
 * 
 * PARTE COMUN DE LA PRACTICA
 */

#ifndef _COMUN_H
#define _COMUN_H

/*INCLUDES*/

#include <vector>
#include <iostream>

/////////////////////END INCLUDES////////////////////


/*DEFINES, CONSTANTS AND GLOBAL FUNCTIONS*/

//#define Solucion pair<vector<int>, unsigned long> 

using namespace std;

const int MAX_EVALS = 50000;
const int N_EXECS = 1;
const int MAX_ITERS = 25;

const unsigned long MAX_INT=0xffffffff;
const long long MAX_LONG=0xffffffffff;

/*STRUCT Measure to print the measurement of either execution
 *
 */
struct Measure {
    long long best, worst;
    double mean_time;
    long long mean_cost;

    Measure() {
        worst = 0;
        best = MAX_LONG;
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

struct Solucion {
    vector<int> first;
    long long second = MAX_LONG;

    bool operator<(const Solucion& g2) const {
        return second < g2.second;
    }

    bool operator==(const Solucion& g2) const {
        return first == g2.first && second == g2.second;
    }
};

struct assign {
    int u, l;
    int cost_add = -1;

    bool operator<(const assign &other) {
        return cost_add < other.cost_add;
    }
};

/*Read matrices from file "filep"
 * Suposse that format of file is:
 *	n -> size of square matrix
 *	m1
 *	m2
 */
bool readMatrices(vector<vector<int> > &f, vector<vector<int> > &d, string filep);

/////////////*END DEFINES, CONSTANTS AND GLOBAL FUNCTIONS*/////////////////////

/*Return cost of solution S*/
int cost(const vector<int> &S, const vector<vector<int> > &f, const vector<vector<int> > &d);


//Reset an fill the vector with sequence from 1 to nelem

void fillVector(vector<int> &v, int tam);

///////////////////////UTILS BL/////////////////////////////////

/*Cost-Fact: cost of swap r and s in S*/
int costFact(const vector<int> &S, const vector<vector<int> > &f, const vector<vector<int> > &d, int r, int s);

/*Operator neighbor: obtain a new neighbor*/
void Neigh(vector<int> &S, int i, int j);

/*Local Search*/
Solucion BL(const Solucion &ini, const vector<vector<int> > &f, const vector<vector<int> > &d);

/*Greedy*/
Solucion Greedy(const vector<vector<int> > &f, const vector<vector<int> > &d);

////////////////////*END UTILS BL*//////////////////////////////

//////////////////////////*ALGORITHMS*////////////////////////////

Solucion ES(const Solucion &s, vector<vector<int> > f, vector<vector<int> > d);

Solucion GRASP(const vector<vector<int> > &f, const vector<vector<int> > &d);

Solucion ILS(const vector<vector<int> > &f, const vector<vector<int> > &d);

Solucion BMB(const vector<vector<int> > &f, const vector<vector<int> > &d);

Solucion ILS_ES(const vector<vector<int> > &f, const vector<vector<int> > &d);
///////////////////////*END ALGORITHMS*///////////////////////////

#endif


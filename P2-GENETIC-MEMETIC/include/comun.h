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

//#define genoma pair<vector<int>, unsigned long> 

using namespace std;

const int MAX_EVALS = 50000;
const int N_EXECS = 1;
const int SIZE_POP = 50;

const double PROB_CROSS_AGG = 0.7;
const double PROB_CROSS_AGE = 1;
const double PROB_MUT = 0.001;
const int K_TOUR = 2;

const int MAX_EVAL_NEIGH = 400;
const double PROB_POP_MEM = 0.1;
const int GENERATIONS_TO_MEM = 10;

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
bool readMatrices(vector<vector<int> > &f, vector<vector<int> > &d, string filep);

/////////////*END DEFINES, CONSTANTS AND GLOBAL FUNCTIONS*/////////////////////

/*Return cost of solution S*/
int cost(const vector<int> &S, const vector<vector<int> > &f, const vector<vector<int> > &d);


//Reset an fill the vector with sequence from 1 to nelem

void fillVector(vector<int> &v, int tam) ;

/*GENETICS ALGORITHMS*/

struct genoma {
	vector<int> first;
	unsigned long second = 0;

	bool operator<(const genoma& g2) const {
		return second < g2.second;
	}

	bool operator==(const genoma& g2) const {
		return first == g2.first && second == g2.second;
	}
};

///////////////////////UTILS BL/////////////////////////////////

/*Cost-Fact: cost of swap r and s in S*/
int costFact(const vector<int> &S, const vector<vector<int> > &f, const vector<vector<int> > &d, int r, int s) ;

/*Operator neighbor: obtain a new neighbor*/
void Neigh(vector<int> &S, int i, int j) ;

/*Local Search*/
genoma BL(const vector<int> &ini, const vector<vector<int> > &f, const vector<vector<int> > &d) ;

////////////////////*END UTILS BL*//////////////////////////////


//other functions
//get best crom

void getBest(const vector<genoma>& pop, genoma &best_gen, int &best) ;

//selection by tournament

int tournament(int k, const vector<genoma>& pop_curr, vector<int>& parents) ;

/* preserve elitism: replace the best by the worst in new population 
 * precondition: best genoma before calculate
 */
void preserveElitism(const genoma& best, vector<genoma>& new_pop) ;

//end other functions

//evaluate pop

void evalPop(vector<genoma>& pop, const vector<vector<int> > &f, const vector<vector<int> > &d) ;

//Start population

void generateStartPop(int size_gen, vector<genoma>& pop, int &best, const vector<vector<int> > &f, const vector<vector<int> > &d) ;


//Selectg

void selectg(const vector<genoma>& pop_curr, vector<int>& parents) ;

//Selects

void selects(const vector<genoma>& pop_curr, vector<int>& parents) ;

//Operators

/* Mut operator: swap gen with other random pos */
void mutOp(genoma& S, const vector<vector<int> > &f, const vector<vector<int> > &d) ;

/* Cross operator: based-position
 * return true if best is replaced, else false
 */
bool crossOpPos(const genoma &par1, const genoma &par2, vector<genoma>& offsprings) ;

/* Cross operator:PMX
 * return true if best is replaced, else false
 */
bool crossOpPMX(const genoma &par1, const genoma &par2, vector<genoma>& offsprings) ;

/*run genetic generational algorithm*/
genoma AGG(const vector<vector<int> > &f, const vector<vector<int> > &d, int size_gen, int TYPE_AGG) ;

/*stationary genetic*/
genoma AGS(const vector<vector<int> > &f, const vector<vector<int> > &d, int size_gen, int TYPE_AGG) ;

/*run memetic algorithm*/
genoma AM(const vector<vector<int> > &f, const vector<vector<int> > &d, int size_gen, int TYPE_AGG, int TYPE_AM) ;

///////////////////////*END GENETIC ALGORITHMS*///////////////////////////

#endif


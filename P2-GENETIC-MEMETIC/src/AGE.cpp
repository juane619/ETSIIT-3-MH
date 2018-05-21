/* 
 * File:   AGG.cpp
 * Author: juane
 *
 * Created on 07 March 2018, 23:44
 * AGG
 * USO: ./agg [<"name.dat">] <seed> <1,2>
 */

#include <dirent.h>
#include <fstream>
#include <algorithm>
#include "comun.h"
#include "timer.h"
#include "random.h"

int main(int argc, char** argv) {
	if (argc > 4 || argc < 3) { // no run
		cerr << "Uso: ./bl [<file.data>] <seed> <TYPE_AG>" << endl;
		return -1;
	} else if (argc == 3 || argc == 4) {
		cout << "RUNNING AG GENERACIONAL ";
		int tam = 0;
		vector<vector<int> > f;
		vector<vector<int> > d;
		double time_elapsed = -1;

		Measure med_bl = Measure();
		ofstream filout("meds");

		if (argc == 3) { // run pasando solo semilla
			unsigned long seed = stoul(argv[1]);
			int TYPE_AGG = atoi(argv[2]);
			
			if(TYPE_AGG==1){
				cout << " OPERADOR DE POSICION\n\n";
			}else if(TYPE_AGG==2){
				cout << " OPERADOR PMX\n\n";
			}

			// obtain a time-based seed:
			//unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();	//Seed establecida a traves de la hora actual

			//GUARDAMOS ARCHIVOS .dat
			vector<string> archivos;
			DIR *dir;
			dirent *ent;

			if (DIR * dir = opendir("instancias/")) {
				while (dirent * ent = readdir(dir)) {
					string archivo(ent->d_name);
					if (archivo.rfind(".dat") != std::string::npos) {
						string nameout("instancias/");
						nameout += archivo;
						archivos.push_back(nameout);
					}
				}
				closedir(dir);
			}

			//

			//Para cada archivo .dat, ejecutamos programa
			sort(archivos.begin(), archivos.end());

			for (auto elem : archivos) {
				Set_random(seed);

				if (readMatrices(f, d, elem)) {
					//output file to boxplot
					//					string name(elem);
					//
					//					string nameout("meds/");
					//					nameout += name;

					//ofstream filout(nameout);

					tam = f.size();

					//EXEC ALGORITHM

					for (int i = 0; i < N_EXECS; i++) {
						start_timers();
						genoma best = AGG(f, d, tam, TYPE_AGG);
						time_elapsed = elapsed_time();

						//						if (best.second > med_bl.worst)
						//							med_bl.worst = best.second;
						//						if (best.second < med_bl.best)
						//							med_bl.best = best.second;
						med_bl.mean_time += time_elapsed;
						med_bl.mean_cost += best.second;

						//cout << sol_BL.second << endl;
					}

					if (N_EXECS > 1) {
						med_bl.mean_time = med_bl.mean_time / N_EXECS;
						med_bl.mean_cost = med_bl.mean_cost / N_EXECS;

						med_bl.print();
					} else {
						filout << med_bl.mean_cost << " " << med_bl.mean_time << endl;
						cout << "INSTANCIA: " << elem << "; mejor caso: " << med_bl.mean_cost << "; tiempo: " << med_bl.mean_time << endl;
						//cin.get();
					}

					cout << endl;

					med_bl = Measure();
					time_elapsed = -1;
				} else
					cerr << "Can't open the file " << elem << "\n";
			}

		} else if (argc == 4) { // run pasando archivo y semilla
			unsigned long seed = stoul(argv[2]);
			int TYPE_AGG = atoi(argv[3]);

			// obtain a time-based seed:
			//unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();	//Seed establecida a traves de la hora actual

			Set_random(seed);

			//Para cada archivo .dat, ejecutamos programa
			if (readMatrices(f, d, argv[1])) {
				//output file to boxplot
				string name(argv[1]);
				name = name.substr(11);
				string nameout("meds/");
				nameout += name;

				tam = f.size();

				for (int i = 0; i < N_EXECS; i++) {
					start_timers();
					genoma best = AGG(f, d, tam, TYPE_AGG);
					time_elapsed = elapsed_time();

					//					if (best.second > med_bl.worst)
					//						med_bl.worst = best.second;
					//					if (best.second < med_bl.best)
					//						med_bl.best = best.second;
					med_bl.mean_time += time_elapsed;
					med_bl.mean_cost += best.second;

					//cout << sol_BL.second << endl;
					//filout << sol_BL.second << endl;
				}
				//filout.close();

				if (N_EXECS > 1) {
					med_bl.mean_time = med_bl.mean_time / N_EXECS;
					med_bl.mean_cost = med_bl.mean_cost / N_EXECS;


					med_bl.print();
				} else {
					cout << "INSTANCIA: " << name << "; mejor caso: " << med_bl.mean_cost << "; tiempo: " << med_bl.mean_time << endl;
				}

				cout << endl;
				//
				//        for (int i = 0; i < sol_BL.first.size(); i++) {
				//            cout << sol_BL.first[i] << " ";
				//        }
				//	cout << endl << sol_BL.second << endl;

			} else
				cerr << "Can't open the file " << argv[1] << "\n";
		}

		filout.close();
	}

	return 0;
}



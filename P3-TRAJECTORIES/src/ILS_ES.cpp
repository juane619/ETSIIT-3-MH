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

#include "timer.h"
#include "random.h"
#include "comun.h"

int main(int argc, char** argv) {
	if (argc < 2 || argc > 3) { // no run
		cerr << "Uso: ./qap <seed> [<file.data>]" << endl;
		return -1;
	} else if (argc == 2 || argc == 3) {
		vector<vector<int> > f;
		vector<vector<int> > d;

		double time_elapsed = -1;

		Measure med_bl = Measure();
		ofstream filout("meds");

		unsigned long seed = stoul(argv[1]);
		// obtain a time-based seed:
		//unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();	//Seed establecida a traves de la hora actual

		//Switch tipo algoritmo
		cout << "RUNNING ILS_ES..\n";

		if (argc == 2) { // run pasando solo semilla
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
			}else{
				cerr << "No es posible abrir el directorio..\n";
				return -1;
			}

			//

			//Para cada archivo .dat, ejecutamos programa
			sort(archivos.begin(), archivos.end());

			for (auto elem : archivos) {
				Set_random(seed);
				srand(seed);

				if (readMatrices(f, d, elem)) {
					cout << elem << endl;
					//output file to boxplot
					//					string name(elem);
					//
					//					string nameout("meds/");
					//					nameout += name;

					//ofstream filout(nameout);


					//Sol inicial
					Solucion s_ini;
					fillVector(s_ini.first, f.size());

					//EXEC ALGORITHM

					for (int i = 0; i < N_EXECS; i++) {
						start_timers();
						//Solucion best = ES(s_ini, f, d);
						Solucion best = ILS_ES(f, d);
						time_elapsed = elapsed_time();

						med_bl.mean_time += time_elapsed;
						med_bl.mean_cost += best.second;

						//cout << best.second << endl;
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
		} else if (argc == 3) { // run pasando archivo y semilla
			Set_random(seed);
			srand(seed);
			string name(argv[2]);

			//Para cada archivo .dat, ejecutamos programa
			if (readMatrices(f, d, name)) {
				//output file to boxplot
				name = name.substr(11);
				string nameout("meds/");
				nameout += name;

				//Sol inicial
				Solucion s_ini;
				fillVector(s_ini.first, f.size());

				for (int i = 0; i < N_EXECS; i++) {
					start_timers();
					//Solucion best = ES(s_ini, f, d);
					Solucion best = ILS_ES(f, d);
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

			} else
				cerr << "Can't open the file " << name << "\n";
		}

		filout.close();
	}

	return 0;
}



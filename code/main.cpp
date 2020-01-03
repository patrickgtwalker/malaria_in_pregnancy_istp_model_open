#define _CRT_SECURE_NO_WARNINGS
#include <sys/stat.h>
#include <omp.h>

#include "library.h"
#include "library/standard.cpp"
#include "find_params.cpp"

#include "Storage.h"
#include "Gen_params.h"
#include "Simulate_pregnancy.h"
#include "simulation.h"

#include "store_params.cpp"
#include "Gen_parms.cpp"
#include "Simulate_pregnancy.cpp"
#include "simulation.cpp"

using namespace std;
#define HOME

#ifdef HOME
int main(int argc, char *argv[]) {
	const clock_t begin = clock();
	const string root = "C:/Users/pgw06/source/repos/MiP_model/Release/";
//const string root = "C:/Users/pgw06/Documents/Model_with_linked_RDT/required/";
	const string directory = root + "def_direct.txt";
	const string name = root + "output/check_recalc2.txt";

#else
int main(int argc, char *argv[]) {
	/// initialise clock for timing
	const clock_t begin = clock();
	/* check there are the minimum files:
	1) root to the model files
	2) directory of lists of model parameters
	3) name for the output
	*/
	int min_arg = 3;
	if (argc - 1 < min_arg)
		error_crit("There should be at least " + as_string(min_arg) + " arguments on the command line, there were " + as_string(argc - 1));
	int arg = 1;
	const string root = as_string(argv[arg++]) + "/";
	const string directory = root + as_string(argv[arg++]);
	const string name = root + as_string(argv[arg++]);
#endif
	// open a cout file to summarise simulation and parameters
	
	string cout_name = name;
	cout_name.replace(cout_name.find(".txt"), 4, ".cout");
	ofstream out_cout;
	out_cout.open(cout_name.c_str(), ios::out);
	cout.rdbuf(out_cout.rdbuf());
	
	/* get the names of the files with for
	1) simulation options
	2) pregnancy_model_parameters
	3) fertility rates
	*/
	map<string, string> file_map;
	ifstream input(directory.c_str());
	if (!input)
		error_crit("Can't open file " + directory + " for input");
	while (!input.eof()) {
		string s;
		string t;
		input >> s >> t;
		if (!input.eof()) {
			if (!input.good())
				error_crit("Error in importing data from " + directory);
			file_map.insert(pair<string, string>(s, t));
		}
	}
	input.close();
	const string sim_params = root + file_map["sim_params"];
	const string preg_params = root + file_map["preg_params"];
	const string fertility_rates = root + file_map["fertility_rates"];
#ifndef HOME
	// log name of executable
	cout << "executable\t" << argv[0] << '\n';
	// read in any parameters changes on the command line
	int arg1 = arg;
	while (arg1 < argc) {
		const string s = as_string(argv[arg1++]);
		if (s == "add")
			continue;
		if (arg1 == argc)
			error_crit(
				"Not counting the word add, there should be an even number of additional command line arguments, name1 value1 name2 value2 ...,  after the first "
				+ as_string(min_arg));
		const string f = as_string(argv[arg1++]);
		const int i = s.find("file");
		if (i != string::npos) {
			if (file_map.find(s) == file_map.end())
				error_crit(s + " is not one of the input list of files");
			file_map[s] = f;
		}
	}
#endif
	// cout date/root/directory used
	cout << '\n' << datestr() << '\t' << timestr() << '\n';
	cout << "executable\t" << argv[0] << '\n';
	cout << "root\t" << root << '\n';
	cout << "directory\t" << directory << '\n';
	cout << "\noutput file\t" << name << "\n";
	////IMPORT SIMULATION PARAMETERS///
	cout << "\nsim parameters\t" << sim_params << '\n';
	cout << "\npreg parameters\t" << preg_params << '\n';
	import_to_map(sim_params, true);
	import_to_map(preg_params, true);
#ifndef HOME 
	if (arg < argc) {
		cout << '\n';
		cout << "\nParameters modified on the command line\n";
	}
	bool to_add = false;
	while (arg < argc) {
		const string s = as_string(argv[arg++]);
		if (s == "add") {
			to_add = true;
			continue;
		}
		if (arg == argc)
			error_crit(
				"Not counting the word add, there should be an even number of additional command line arguments, name1 value1 name2 value2 ...,  after the first "
				+ as_string(min_arg));
		const string f = as_string(argv[arg++]);
		const int i = s.find("file");
		if (i == string::npos) {
			const double x = atof(f.c_str());
			to_add ? add_to_map(s, x) : modify_map(s, x);
			cout << s << '\t' << x << endl;
		}
	}
#endif
	////IMPORT FERILITY RATES ///
	cout << "\nfertility_rates\t" << fertility_rates << '\n';
	vector<vector<double>> f_rates = store_rates(fertility_rates);
	simulation simulation;
	simulation.ratesarray=store_rates(fertility_rates);
	simulation.summary = from_map_bool("summary", 0);
	simulation.file.open(name);
	if (simulation.summary) {
		string summary_name = name;
		summary_name.replace(summary_name.find(".txt"), 4, "_summary.txt");
		simulation.file_summary.open(summary_name);
	}
	simulation.num_sims = from_map("num_women", 1, 10000000000000000000, -99999);
	////IMPORT AND SET UP OUTPUT CATEGORIES///
	int output_strats = 0;
	vector<double> output_cats_0;
	vector<double> output_cats_1;
	while (in_map("output_beg_" + as_string(output_strats))) {
		double par_beg = from_map("output_beg_" + as_string(output_strats), -1, 10000, 1);
		double par_end = from_map("output_end_" + as_string(output_strats), -1, 10000, 1);
		simulation.par_down.push_back(par_beg >= 0 ? par_beg : 0);
		simulation.par_up.push_back(par_end >= 0 ? par_end : 200);
		output_strats++;
	}
	////SET UP DETERMINISTIC GENERAL POP MODEL///
	simulation.gen_parms.setup();
	simulation.gen_parms.ft = from_map("ft", 0, 1, 1);
	simulation.gen_parms.EIR = from_map("EIR", 0, 10000, 1);
	simulation.gen_parms.init_determ_model();
	////RUN SIMULATED PREGNANCY COHORT///
	simulation.run_simulation();
	cout << (double)(clock() - begin) / CLOCKS_PER_SEC << " seconds \n";
	return 0;
}
		
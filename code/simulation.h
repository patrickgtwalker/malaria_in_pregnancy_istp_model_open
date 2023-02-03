/// FOR SETTING UP, RUNNING THE PREGNANCY COHORT AND STORE OUTPUT


struct simulation{
	int num_sims;
	bool summary;
	bool inf_history;
	bool EIR_vector;
	gen_parms gen_parms;
	pregnancy preg;
	vector<double> EIRs_to_run;
	vector<vector<double>> ratesarray;
	vector<double> par_down;
	vector<double> par_up;
	vector<double> num_simulated;
	vector<vector<double>> week_prof;
	vector<vector<double>> week_peri_prof;
	vector<vector<double>> previous_inf_dist;
	vector<vector<double>> previous_inf_dist_inf;
	double primi_prev_anc1;
	ofstream file;
	vector<double> prop_peri;
	vector<double> prop_plac;
	vector<double> peri_dur;
	vector<double> plac_dur;
	vector<double> fail_cure;
	vector<double> prop_fail_proph;
	vector<double> tot_fail_proph;
	vector<double> LBW_risk;
	ofstream file_summary;
	ofstream file_inf_history;


	void run_simulation(void);
	void setup_summary(void);
	void setup_inf_history(void);
	void setup_timeline(void);
	void ANC_setup(void);
	
	void life_time(void);
	void run_pregnancy(pregnancy& run_preg);
	void getfirststate(void);

	void store_output(pregnancy& store_preg);
	void store_summary(pregnancy& store_preg);
	void store_hist_inf_dist(pregnancy& store_preg);
	void write(void);
	void write_summary(void);
	void write_inf_history(void);
	void clear_sim(void);
};

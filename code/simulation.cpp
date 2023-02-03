


void simulation::run_simulation(void) {
	if (!EIR_vector) EIRs_to_run.push_back(gen_parms.EIR);
	if (summary)setup_summary();
	if (inf_history)setup_inf_history();
	for (int i = 0;i < EIRs_to_run.size();i++) {
		cout << EIRs_to_run[i] << "\n";
		clear_sim();
		gen_parms.EIR = EIRs_to_run[i];
		gen_parms.init_determ_model();
		setup_timeline();
		ANC_setup();
		for (int j = 0; j < num_sims; j++) {
			preg.clearall();
			life_time();
		}
		//write();
		//if (summary)write_summary();
		if (inf_history)write_inf_history();
	}
	return;
}


void simulation::setup_timeline(void){
	week_prof.resize(par_down.size(),vector<double>(280,0));
	week_peri_prof.resize(par_down.size(),vector<double>(280,0));
	num_simulated.resize(par_down.size());
	return;
}

void simulation::setup_summary(void) {
	prop_peri.resize(par_down.size(),0);
	prop_plac.resize(par_down.size(), 0);
	peri_dur.resize(par_down.size(), 0);
	plac_dur.resize(par_down.size(), 0);
	fail_cure.resize(par_down.size(), 0);
	prop_fail_proph.resize(par_down.size(), 0);
	tot_fail_proph.resize(par_down.size(), 0);
	if (preg.strategy == 0) LBW_risk.resize(par_down.size(), 0);
	return;
}

void simulation::setup_inf_history(void) {
	primi_prev_anc1 = 0;
	previous_inf_dist.resize(par_down.size(),vector<double> (20,0));
	previous_inf_dist_inf.resize(par_down.size(), vector<double>(20, 0));
	file_inf_history << "EIR\tprimi_ANC_prev";
	for (int i = 0; i < par_down.size(); i++) {
		for (int j = 0; j < previous_inf_dist_inf[i].size(); j++) {
			file_inf_history << "\tprev_inf_dist_inf_cat_" + as_string(i) + "_n_" + as_string(j);
		}
		for (int j = 0; j < previous_inf_dist[i].size(); j++) {
			file_inf_history << "\tprev_inf_dist_cat_" + as_string(i) + "_n_" + as_string(j);
		}
	}
	file_inf_history << "\n";
	return;
}

void simulation::ANC_setup(void) {
	preg.ft = gen_parms.ft;
	preg.strategy = from_map("strategy", 0, 4, 1);
	preg.intervene = from_map("intervene", 0, 1, 1);
	preg.IPTeff = from_map("drug_efficacy_" + as_string(from_map("IPT_drug", 1, 2)), 0, 1);
	double IPThl = from_map("drug_half_life_" + as_string(from_map("IPT_drug", 1, 2)), 0, 10000);
	preg.IPTshape = from_map("drug_shape_" + as_string(from_map("IPT_drug", 1, 2)), 0, 10000);
	preg.IPTscale = pow(IPThl / tgamma(1 + 1 / preg.IPTshape), -preg.IPTshape);
	preg.ISTeff = from_map("drug_efficacy_" + as_string(from_map("IST_drug", 1, 2)), 0, 1);
	double ISThl = from_map("drug_half_life_" + as_string(from_map("IST_drug", 1, 2)), 0, 10000);
	preg.ISTshape = from_map("drug_shape_" + as_string(from_map("IST_drug", 1, 2)), 0, 10000);
	preg.ISTscale = pow(ISThl / tgamma(1 + 1 / preg.ISTshape), -preg.ISTshape);
	preg.failprob = 1 - parameter_map["efficacy_previous_failure_deduction"].first;
	preg.perfect_test = from_map_bool("perfect_test", 0);
	preg.perfect_sens = from_map("perfect_sens", 0, 1);
	preg.prog1 = from_map("prog_acute", 0, 100000);
	preg.prog2 = from_map("prog_chron", 0, 100000);
	preg.Kp = from_map("sequest_imm_power", 0, 100000);
	preg.Ip0 = from_map("sequest_imm_offset", 0, 100000);
	preg.Kchron = from_map("chron_imm_power", 0, 100000);
	preg.Ichron0 = from_map("chron_imm_offset", 0, 100000);
	preg.rate = from_map("LBW_hazard", 0, 100000);
	int ANC_visits = 1;
	while (in_map("ANC_" + as_string(ANC_visits)) && from_map("ANC_" + as_string(ANC_visits)) == 1) {
		preg.ANC_times.push_back(from_map("ANC_time_" + as_string(ANC_visits), 0, 280));
		ANC_visits++;
	}
	preg.ru = 1.00 / 107.867;
	preg.rd = 1.00 / 5.00;
	preg.rp = 1 / 9.00;
	preg.inflag = 7.0;
	preg.gest_time = 0;
	preg.gestation_duration = 40.0 * 7.0;
	preg.blood_to_IVS = 12.0 * 7.0;
	preg.histplac = 0;
	preg.parity = 0;

	for (int k = 0; k < 200; k++) {
		preg.progC[k] = preg.prog2 * (1 + pow(k / preg.Ichron0, preg.Kchron));
		preg.fracA[k] = 1 / ((1 + pow(k / preg.Ip0, preg.Kp)));
	}


	preg.deathprob = 0.009208788 / 365.0;
	preg.weekinf.resize(preg.gestation_duration, 0);
	preg.weekperiinf.resize(preg.gestation_duration, 0);
	preg.primi_sens = from_map("primi_sens", 0, 10000);
	preg.rdt_preg_shape = from_map("rdt_preg_shape", 0, 10000);
	preg.rdt_preg_offset = from_map("rdt_preg_offset", 0, 10000);
	preg.rdt_adult_shape = from_map("rdt_adult_shape", -10000, 10000);
	preg.rdt_adult_offset = from_map("rdt_adult_offset", -10000, 10000);
	preg.intervene_imm = from_map_bool("intervene_imm", 1);
	preg.late_PG_sens_odds = from_map("late_PG_sens_odds", -10000, 10000);
	preg.late_grav_OR = from_map("late_grav_OR", -10000, 10000);
	preg.prev_inf_OR = from_map("prev_inf_OR", -10000, 10000);
	preg.first_tri_rdt = from_map("first_tri_rdt", 0, 2);
	preg.first_tri_visit = from_map_bool("first_tri_visit", 0);
	preg.first_tri_visit_time = from_map("first_tri_visit_time", -1, 84);
	double num_sims = parameter_map["num_women"].first;
	preg.sensitivity.resize(20, 0.0);
	double non_preg_rdt = 1 / (1 + exp(-(preg.rdt_adult_shape * log(gen_parms.pcr_adult / (1 - gen_parms.pcr_adult)) + preg.rdt_adult_offset)));
	for (int par = 0; par < preg.sensitivity.size(); par++) {
		double preg_odds = 1 / ((gen_parms.pcr_adult/ non_preg_rdt) - 1) * (1 + preg.primi_sens / (1 + pow(par / preg.rdt_preg_offset, preg.rdt_preg_shape)));
		preg.sensitivity[par] = preg_odds / (1 + preg_odds);
	}
	preg.npreg_rdt_sens = non_preg_rdt / gen_parms.pcr_adult;
	for (int i = 0; i < 20; i++) {
		preg.late_sensitivity.push_back(exp(preg.late_PG_sens_odds + preg.late_grav_OR * i) / (1 + exp(preg.late_PG_sens_odds + preg.late_grav_OR * i)));
		preg.late_sensitivity_prev_inf.push_back(exp(preg.late_PG_sens_odds + preg.late_grav_OR * i + preg.prev_inf_OR) / (1 + exp(preg.late_PG_sens_odds + preg.late_grav_OR * i + preg.prev_inf_OR)));
	}
	return;
}

void simulation::life_time(void) {
		//reset lifetime
		preg.parity = 0;
		preg.histplac = 0;
		preg.histinf = 0;
		//find heterogeneity category
		double rand = runif();
		double hetsum = 0;
		int hetcount = -1;
		while (hetsum < rand) {
			hetcount++;
			hetsum += gen_parms.het_wt[hetcount];
		}
		preg.hetcat = hetcount;
		//simulate lifetime for age_cat 25 through 32 (15-49), account for natural mortality
		double trackdate = 0;
		double death = 0;
		for (int j = 25; j < 32; j++) {
			preg.age_cat = j;
			double next_age_cat = 5.0 * (j - 24.0) * 365.0;
			while (trackdate < next_age_cat && death == 0) {
				double birth_time = genexp(ratesarray[j - 25][preg.parity]) * 365.0;
				if ((birth_time + trackdate) < next_age_cat) {
					double prob = runif();
					if (prob > (1 - exp(-(birth_time + preg.gestation_duration) * preg.deathprob))) {
						/// run the pregnancy
						run_pregnancy(preg);
						trackdate = trackdate + birth_time + preg.gestation_duration;
						preg.parity++;
					}
					else { death++; }
				}
				else {
					trackdate = next_age_cat;
				}
			}
		}
		return;
	}




void simulation::run_pregnancy(pregnancy& run_preg) {
	//// RESET FROM PREVIOUS PREGNANCY ////
	run_preg.clearall();
	//// FIND INFECTION STATE AT CONCEPTION ///
	//simulation.get_first_state(this_preg);
	getfirststate();
	/// GENERATE TIMES OF EXPOSURE ////
	run_preg.generateperiinfs();
	/* intervene ==1 is to look at what happens when an intervention is
	newly implemented so women have not yet received IPTp in a previous pregnancy
	to do this we simulate each pregnancy twice, with the same underlying exposure,
	one to look the impact of the intervention in a women with gravidity x,
	one to generate immunity for a woman with gravidity x+1 who hasn't received IPTp before
	*/
	if ((!run_preg.intervene_imm) & (run_preg.strategy != 0)) {
		pregnancy intervene_preg = run_preg;
		intervene_preg.run_to_delivery();
		// store the pregnancy with the intervention
		store_output(intervene_preg);
		if(summary)store_summary(intervene_preg);
		if (inf_history)store_hist_inf_dist(intervene_preg);
		int this_strat = run_preg.strategy;
		//run the pregnancy without the intervention
		run_preg.strategy = 0;
		run_preg.run_to_delivery();
		run_preg.strategy = this_strat;
	}
	else {
		// no need to run two pregnancies 
		run_preg.run_to_delivery();
		store_output(run_preg);
		if (summary)store_summary(run_preg);
		if (inf_history)store_hist_inf_dist(run_preg);
	}
	return;
}


/// FUTURE PERIPHERAL DYNAMICS DUE TO INFECTION STATUS AT CONCEPTION
void simulation::getfirststate(void) {
	double strand = runif();
	preg.pcr_prev = (gen_parms.D[preg.age_cat][preg.hetcat] + gen_parms.A[preg.age_cat][preg.hetcat] + gen_parms.U[preg.age_cat][preg.hetcat]) / gen_parms.den_het[preg.age_cat][preg.hetcat];
	preg.recA = gen_parms.recI[preg.age_cat];
	preg.cFOI = gen_parms.FOI[preg.age_cat][preg.hetcat];
	preg.fracD = gen_parms.phis[preg.age_cat][preg.hetcat];
	//cout << "check gen " << gen_parms.FOI[25][2] << "\n";
	//cout << " rp " << preg.rp << " rd " << preg.rd << " ru " << preg.ru << " rp " << preg.rp << " recA " << preg.recA << " age cat "<< preg.age_cat<<"\n";
	preg.init_plac = 0;
	if (gen_parms.T[preg.age_cat][preg.hetcat] / gen_parms.den_het[preg.age_cat][preg.hetcat] > strand) {
		/// CURRENTLY TAKING TREATMENT - incorporate future prophylaxis
		preg.gest_time = genexp(preg.rd) + genexp(preg.rp);
	}
	else if ((gen_parms.T[preg.age_cat][preg.hetcat] + gen_parms.P[preg.age_cat][preg.hetcat]) / gen_parms.den_het[preg.age_cat][preg.hetcat] > strand) {
		/// CURRENTLY PROTECTED DUE TO POST-TREATMENT PROPHYLAXIS- incorporate current prophylaxis
		preg.gest_time = genexp(preg.rp);
	}
	else if ((gen_parms.T[preg.age_cat][preg.hetcat] + gen_parms.P[preg.age_cat][preg.hetcat] + gen_parms.D[preg.age_cat][preg.hetcat] + gen_parms.A[preg.age_cat][preg.hetcat] + gen_parms.U[preg.age_cat][preg.hetcat])/ gen_parms.den_het[preg.age_cat][preg.hetcat] > strand) {
		//CURRENTLY INFECTED
		preg.ever_peri = 1;
		double cleartime;
		if ((gen_parms.T[preg.age_cat][preg.hetcat] + gen_parms.P[preg.age_cat][preg.hetcat] + gen_parms.D[preg.age_cat][preg.hetcat]) / gen_parms.den_het[preg.age_cat][preg.hetcat] > strand) cleartime = genexp(preg.rd) + genexp(preg.recA) + genexp(preg.ru);/* CURRENTLY WITH CLINICAL DISEASE THAT WILL NOT BE TREATED */
		else if ((gen_parms.T[preg.age_cat][preg.hetcat] + gen_parms.P[preg.age_cat][preg.hetcat] + gen_parms.D[preg.age_cat][preg.hetcat] + gen_parms.A[preg.age_cat][preg.hetcat]) / gen_parms.den_het[preg.age_cat][preg.hetcat] > strand) cleartime = genexp(preg.recA) + genexp(preg.ru);/* CURRENTLY WITH RECENT HIGH-DENSITY ASYMPTOMATIC - clear sooner*/
		else cleartime = genexp(preg.ru);/*CURRENTLY WITH OLDER LOW-DENSITY ASYMPTOMATIC - clear soonest*/
		// STORE START AND END OF PERIPHERAL INFECTION
		preg.start_peri_inftimes_B.push_back(0);
		preg.end_peri_time_K.push_back(cleartime);
		///SEQUESTERS IF PERSITS UNTIL BLOOD ENTERS IVS AND TIME TO SWITCH TO PLACENTAL PHENOTYPE
		if (cleartime > preg.blood_to_IVS + preg.inflag) {
			preg.init_plac = 1;

		}
		else preg.init_plac = 0;
	}
	return;
}

void simulation::store_hist_inf_dist(pregnancy& store_preg) {
	for (int i = 0; i < par_down.size(); i++) {
		if ((store_preg.parity >= par_down[i]) & (store_preg.parity <= par_up[i])) {
			if (previous_inf_dist[i].size() <=store_preg.hist_inf_beg) {
				previous_inf_dist[i].push_back(1.00);
			}
			else {
				previous_inf_dist[i][store_preg.hist_inf_beg]++;
			}
			if (store_preg.inf_anc1) {
				if (previous_inf_dist_inf[i].size() <=store_preg.hist_inf_beg) {
					previous_inf_dist_inf[i].push_back(1.00);
				}
				else {
					previous_inf_dist_inf[i][store_preg.hist_inf_beg]++;
				}
			}
		}
	}
	if (store_preg.parity == 0 & store_preg.inf_anc1) primi_prev_anc1++;
	return;
}

void simulation::clear_sim(void) {
	for (auto& i : week_prof)
		std::fill(i.begin(), i.end(), 0);
	for (auto& i : week_peri_prof)
		std::fill(i.begin(), i.end(), 0);
	for (auto& i : previous_inf_dist)
		std::fill(i.begin(), i.end(), 0);
	for (auto& i : previous_inf_dist_inf)
		std::fill(i.begin(), i.end(), 0);
	fill(prop_peri.begin(), prop_peri.end(), 0);
	fill(prop_plac.begin(), prop_plac.end(), 0);
	fill(plac_dur.begin(), plac_dur.end(), 0);
	fill(fail_cure.begin(), fail_cure.end(), 0);
	fill(prop_fail_proph.begin(), prop_fail_proph.end(), 0);
	fill(tot_fail_proph.begin(), tot_fail_proph.end(), 0);
	fill(num_simulated.begin(), num_simulated.end(), 0);
	primi_prev_anc1 = 0;
	return;
}

void simulation::store_output(pregnancy& store_preg) {
	for (int i = 0; i <par_down.size(); i++) {
		if ((store_preg.parity >= par_down[i]) & (store_preg.parity <= par_up[i])) {
			num_simulated[i]++;
			for (int k = 0; k < store_preg.weekinf.size(); k++) {
				week_prof[i][k] += store_preg.weekinf[k];
				week_peri_prof[i][k] += store_preg.weekperiinf[k];

			}
		}
	}
	return;
}
void simulation::store_summary(pregnancy& store_preg) {
	for (int i = 0; i < par_down.size(); i++) {
		if ((store_preg.parity >= par_down[i]) & (store_preg.parity <= par_up[i])) {
			prop_peri[i] += store_preg.ever_peri;
			prop_plac[i] += store_preg.ever_plac;
			plac_dur[i] += store_preg.timeinfected;
			fail_cure[i] += store_preg.clear_fail;
			prop_fail_proph[i] += store_preg.proph_fail;
			tot_fail_proph[i] += store_preg.proph_fail_tot;
			if (preg.strategy == 0) LBW_risk[i]+=1-exp(-store_preg.timechron* store_preg.rate);
		}
	}
	return;
}

void simulation::write_summary(void) {
	file_summary << "Grav_cat\tprop_peri\tprop_plac\taverage_plac_duration\tprop_fail_clear\tprop_newly_infected\taverage_new_infections\t";
	if (preg.strategy == 0)file_summary << "LBW_risk\n"; else file_summary << "\n";
	for (int j = 0; j < par_down.size(); j++) {
		if (par_down[j] == 0 && par_up[j] == 200) file_summary << "All\t";
		else if (par_up[j] == 200) file_summary << as_string(par_down[j]) + "_max\t";
		else file_summary << as_string(par_down[j]) + "_" + as_string(par_up[j]) + "\t";
		file_summary << prop_peri[j] / num_simulated[j] << "\t" << prop_plac[j] / num_simulated[j] << "\t" << plac_dur[j] / num_simulated[j] << "\t" << fail_cure[j] / num_simulated[j] << "\t" << prop_fail_proph[j] / num_simulated[j] << "\t" << tot_fail_proph[j] / num_simulated[j] << "\t";
		if (preg.strategy == 0)file_summary << LBW_risk[j] / num_simulated[j] << "\n"; else file_summary << "\n";
	}
	file_summary.close();
	return;
}

void simulation::write_inf_history(void) {
	file_inf_history <<gen_parms.EIR<<"\t"<<primi_prev_anc1/ num_simulated[0];
	for (int i = 0; i < par_down.size(); i++) {
		for (int j = 0; j < previous_inf_dist_inf[i].size(); j++) {
			file_inf_history << "\t"<<previous_inf_dist_inf[i][j] / num_simulated[i];
		}
		for (int j = 0; j < previous_inf_dist[i].size(); j++) {
			file_inf_history << "\t" << previous_inf_dist[i][j] / num_simulated[i];
		}
	}
	file_inf_history << "\n";
	file_inf_history.flush();
	//file_inf_history.close();
	return;
}



void simulation::write(void) {
	file << "Day\t";
	for (int j = 0; j < par_down.size(); j++) {
		if (par_down[j] == 0 && par_up[j] == 200) {
			file << "peri_infected_all\tplac_infected_all\t";
		}
		else if (par_up[j] == 200) {
			file << "peri_infected_" + as_string(par_down[j]) + "_max\t" << "plac_infected_" + as_string(par_down[j]) + "_max\t";
		}
		else {
			file << "peri_infected_" + as_string(par_down[j]) + "_" + as_string(par_up[j]) + "\t" << "plac_infected_" + as_string(par_down[j]) + "_" + as_string(par_up[j]) + "\t";
		}
	}
	for (int j = 0; j < 280; j++) {
		file << "\n";
		file << j << "\t";
		for (int i = 0; i < week_prof.size(); i++) {
			file << week_peri_prof[i][j] / num_simulated[i] << "\t" << week_prof[i][j] / num_simulated[i] << "\t";
		}
	}
	file << "\n";
	file.close();
	return;

}

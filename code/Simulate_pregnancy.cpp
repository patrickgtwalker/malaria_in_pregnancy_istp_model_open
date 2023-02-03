



	
/// RUNS FROM gest_time AT END OF FIRST STATE TO END OF PREGNANCY TO FIND EXPOSURE TO PERIPHERAL INFECTION IN ABSENCE OF INTERVENTION
void pregnancy::generateperiinfs(void){
	//FIND NEXT INFECTION TIME
	gest_time += genexp(cFOI);
	//RUN UNTIL MATERNAL BLOOD ENTERS INTERVILLOUS SPACE (i.e. placenta too young to be infected)
	while (blood_to_IVS > gest_time) {
		double rand = runif(); //random number to determine course of infection
		if (rand < fracD * ft) {//IS SYMPTOMATIC AND SUBSEQUENTLY TREATED 
			for (int i = 0; i < end_peri_time_K.size(); i++) {
				//truncate any earlier, persisting infections
				if (end_peri_time_K[i] > gest_time) {
					end_peri_time_K[i] = gest_time;
				}
			}
			//remove any future placental exposure due to prior infections
			init_plac = 0;
			//move to end of prophylaxis
			gest_time = gest_time + genexp(rd) + genexp(rp);
			
		}
		else {
			//IS AN UNTREATED INFECTION
			double cleartime = 0;
			if (rand < fracD) cleartime = gest_time + genexp(rd) + genexp(recA) + genexp(ru); //IS SYMPTOMATIC AND GOES UNTREATED 	
			else cleartime = gest_time + genexp(recA) + genexp(ru); 	//IS ASYMPTOMATIC
			if (cleartime>blood_to_IVS + inflag) {
				init_plac = 1; ///WILL SEQUESTERS ONCE PLACENTA IS SUSCEPTIBLE (provided no intervening treatment)
			}
			else init_plac = 0; ///CLEARS TOO EARLY TO SEQUESTER
			//store peripheral infection and clearance times
			start_peri_inftimes_B.push_back(gest_time);
			end_peri_time_K.push_back(cleartime);
		}
		///move on to next infection
		gest_time += genexp(cFOI);
	}
	/// any infections present as blood flows to IVS now stored
	if(init_plac==1)plac_exposure_times_P.push_back(blood_to_IVS + inflag);
	/// run for remainder of gestation (as previously but now placenta can be readily infected)
	while (gestation_duration > gest_time) {
		double rand = runif();
		if (rand < fracD * ft) {//IS SYMPTOMATIC AND SUBSEQUENTLY TREATED 
			for (int i = 0; i < end_peri_time_K.size(); i++) {
				if (end_peri_time_K[i] > gest_time) {
					end_peri_time_K[i] = gest_time;
				}
			}
			gest_time = gest_time + genexp(rd) + genexp(rp);
		}
		else {
			double cleartime = 0;
			if (rand < fracD) cleartime = gest_time + genexp(rd) + genexp(recA) + genexp(ru); //IS SYMPTOMATIC AND GOES UNTREATED 	
			else cleartime = gest_time + genexp(recA) + genexp(ru); 	//IS ASYMPTOMATIC
			if (cleartime > gest_time + inflag) {
				plac_exposure_times_P.push_back(gest_time + inflag);
			}
			start_peri_inftimes_B.push_back(gest_time);
			end_peri_time_K.push_back(cleartime);
		}
		gest_time += genexp(cFOI);
	}
	return;
}

/* ONCE WE HAVE PERIPHERAL DYNAMICS IN ABSENCE OF ANC-BASED INTERVENTION FUNCTION INCLUDES PLACENTAL DYNAMICS INCLUDING
ANY ANC-BASED INTERVENTION */
void pregnancy::run_to_delivery(void) {
	//set ANC counters to zero
	bool first_tri_past = false;
	hist_inf_beg = histinf;
	vector<bool> past_ANC(ANC_times.size() , false);
	if (start_peri_inftimes_B.size() == 0)gest_time = gestation_duration + 1.0; //PREGNANCY NOT EXPOSED 
	else {
		// SOMETHING HAPPENS DURING PREGNANCY
		p1 = prog1;
		p2 = progC[histplac];
		nextplacchange = 0;
		nextplactime = 10000;
		// to ensure we reach end of gestation
		start_peri_inftimes_B.push_back(gestation_duration + 1.0);
		end_peri_time_K.push_back(gestation_duration + 1.0);
		plac_exposure_times_P.push_back(gestation_duration + 1.0);
		//find beginning and end of first peripheral infection and remove from vector
		start_pcr_time = start_peri_inftimes_B[0];
		pcr_clear_time = end_peri_time_K[0];
		start_peri_inftimes_B.erase(start_peri_inftimes_B.begin());
		end_peri_time_K.erase(end_peri_time_K.begin());
		// find next event (either peripheral infection clears or placenta is exposed)
		gest_time = min(pcr_clear_time, plac_exposure_times_P[0]);
		if ((ANC_times.size() > 0) && (!past_ANC[0]) && (gest_time > ANC_times[0]) && (start_pcr_time < ANC_times[0]) && (pcr_clear_time > ANC_times[0])) inf_anc1 = true;
		///// CHECK IF WE'VE ACTUALLY PASSED ANY ANC DATES AND IF SO WHETHER WE PROVIDED PROPHYLAXIS
		if(first_tri_visit) any_first_trimester(first_tri_past);
		for (int ANC = 0; ANC < ANC_times.size(); ANC++) {
			IPTISTupdate(ANC, past_ANC);
		}
		//NOW LOOP UNTIL DELIVERY
		while (gest_time < gestation_duration) {
		
			/// OTHERWISE WE HAVE REACHED DELIVERY
			/*while gest_time is <gestation_duration there are two types of events peripheral or placental
			peripheral:
			- A new infection occurs (pcr positivity begins if not ongoing and eventual clearance time is set)
			- An ongoing infection has the potential to sequester in placenta
			- Infection clears (or would in absence of placenta)
			placental:
			- An infection becomes chronic (i.e. accumulation of visible pigment by histology
			- An infection clears (though pigment may remain)
			'week_peri_inf' then tracks peripheral pcr_prevalence (which may be extended by placental infection)
			'week_plac_inf' only tracks if there is an active placental infection
			*/
			if (gest_time == next_peri_time) {
				//If here an event related to peripheral infection
				if (gest_time == plac_exposure_times_P[0]) {
					/// If here it's a peripheral infection potentially sequestering in placenta
					///INDICATOR FOR EVER EXPOSED TO PLACENTAL PARASITES
						/// see if any new infection occurs after first ANC (so as failure of protect)
					if (gest_time > ANC_times[0]) {
						proph_fail = 1;
						proph_fail_tot++;
					}
					/// see if infection sequesters to appreciable levels in placenta
					if (runif() < fracA[histplac]) {
						ever_plac = 1;
						total_plac++;
						newinfpath();///FIND OUT COURSE OF PLACENTAL INFECTION (WHEN IT BECOMES CHRONIC AND WHEN IT CLEARS) IN ABSENCE OF INTERVENTION	
						nextplacchange = nextstate(); /// FIND OUT TYPE OF NEXT PLACENTAL CHANGE
					}
					// NOW REMOVE ANY EXPOSURES WE'VE JUST PASSED
					while (gest_time == plac_exposure_times_P[0]) plac_exposure_times_P.erase(plac_exposure_times_P.begin());
				}
				else if (gest_time == start_peri_inftimes_B[0]) {
					///NEW PERIPHERAL INFECTION - if no ongoing infection then pcr positivity strats
					if (start_pcr_time > start_peri_inftimes_B[0])start_pcr_time = start_peri_inftimes_B[0];
					// UPDATE CLEARANCE TIME
					pcr_clear_time = end_peri_time_K[0];
					// REMOVE THIS PERIPHERAL INFECTION WE'VE JUST PASSED 
					start_peri_inftimes_B.erase(start_peri_inftimes_B.begin());
					end_peri_time_K.erase(end_peri_time_K.begin());
				}
				else if (gest_time == pcr_clear_time) {
					//PCR positivity in absence of placental sequestration comes to an end
					//STORE PCR positivity
					ever_peri = 1;
					fill(weekperiinf.begin() + (int)start_pcr_time, weekperiinf.begin() + (int)pcr_clear_time, 1);
					pcr_clear_time = 281;
					start_pcr_time = 281;
				}
				//FIND NEXT PERIPHERAL EVENTS
				next_peri_time = min({ plac_exposure_times_P[0], start_peri_inftimes_B[0], pcr_clear_time });
			}
			else if (gest_time == nextplactime) {
				/// IF HERE NEXT EVENT IS A STAGE OF AN EXISTING PLACENTAL INFECTION (ever clearing or becoming 'chronic')
				placentalevent();//CARRIES OUT THIS UPDATE 
			}
			/// FIND OUT WHETHER NEXT EVENT IS PERIPHERAL OF PLACENTAL
			gest_time = ((nextplactime > gest_time && nextplactime < next_peri_time) ? nextplactime : next_peri_time);
			///// CHECK IF WE'VE ACTUALLY PASSED ANY ANC DATES AND IF SO WHETHER WE PROVIDED PROPHYLAXIS
			if (first_tri_visit) any_first_trimester(first_tri_past);
			if ((ANC_times.size() > 0) && (!past_ANC[0]) && (gest_time > ANC_times[0]) && (start_pcr_time < ANC_times[0]) && (pcr_clear_time > ANC_times[0])) inf_anc1 = true;
			for (int ANC = 0; ANC < ANC_times.size(); ANC++) {
					IPTISTupdate(ANC, past_ANC);
				}	
		}
		/// STORE EXPOSURE AT THE END OF PREGNANCY
		histplac += total_plac;
		histinf += ever_plac;
		if (totchron > 0) { timechron += 280 - startchrondate; }
		if (totalpara > 0) {
			ever_peri = 1;
			timeinfected += 280 - startinfdate;
			fill(weekperiinf.begin() + (int)startinfdate, weekperiinf.end(), 1);
			fill(weekinf.begin() + (int)startinfdate, weekinf.end(), 1);
		}
		if (start_pcr_time < 280) {
			ever_peri = 1;
			fill(weekperiinf.begin() + (int)start_pcr_time, weekperiinf.end(), 1);
		}
	}
	return;
}

//////////////////////////////////////// PLACENTAL INFECTION IN ABSENCE OF INTERVENTION ///////////////////////////////
/// PLOTS THE COURSE OF ANY NEW PLACENTAL INFECTION
void pregnancy::newinfpath(void){
		// if currently no placental infection set a new start of placental prevalence
	if(totalpara==0){startinfdate=gest_time;}
	// add one to total current infections
	totalpara++;
	// find time infection becomes chronic and save in order of occurence
	double trackcourse=gest_time+genexp(prog1);
	placevec(start_chronic_placenta,trackcourse);
	// find time active infection clears and save in order of occurence
	double durchron=genexp(p2);
	trackcourse+=durchron;
	placevec(end_placenta_active_R,trackcourse);	
	return;
}

/// FUNCTION TO IMPLEMENT PLACENTAL EVENT AND FIND NEXT ONE
void pregnancy::placentalevent(void) {
	if (nextplacchange == 1) {
		//--- PIGMENT APPEARS AS A RESULT OF INFECTION ---////
		if (totchron == 0) {
			//IF NO ONGOING CHRONIC INFECTION BEGINS HERE
			startchrondate = gest_time;
			//ever_chron = 1;
		}
		totchron++;
		//REMOVE THE EVENT
		removetime(start_chronic_placenta);
	}
	else if (nextplacchange == 2) {
		//--- ACTIVE PARASITAEMIA FROM AN INFECTION IS CLEARED---////
		totchron--;
		totalpara--;
		//if end of chronic exposure store duration for link to LBW
		if (totchron == 0) timechron += gest_time - startchrondate; 
		//if all infection cleared store duration and timing of exposure
		if (totalpara == 0) {
			timeinfected += gest_time - startinfdate;
			if(intervene==1){
				fill(weekinf.begin() + (int)startinfdate, weekinf.begin() + (int)gest_time, 1);
				fill(weekperiinf.begin() + (int)startinfdate, weekperiinf.begin() + (int)gest_time, 1);
			}
		}
		//REMOVE THE EVENT
		removetime(end_placenta_active_R);
	}
	/// FIND THE NEXT TYPE OF PLACENTAL INFECTION AND TIME
	nextplacchange = nextstate();
	return;
}

/// FUNCTION TO FIND NEXT PLACENTAL EVENT AND TIME
double pregnancy::nextstate(void) {
	nextplactime = 10000000;
	double nextplacchange = 0;
	if (start_chronic_placenta.size() > 0 && start_chronic_placenta[0] < nextplactime) {
		nextplactime = start_chronic_placenta[0];
		nextplacchange = 1;
	}
	if (end_placenta_active_R.size() > 0 && end_placenta_active_R[0] < nextplactime) {
		nextplactime = end_placenta_active_R[0];
		nextplacchange = 2;
	}
	return nextplacchange;
}


/// FUNCTION TO ADD A NEW TIME IN SEQUENCE
void pregnancy::placevec(vector<double>& vec, double ins_time) {
	int vsize = vec.size();
	if (vsize < 1) { vec.push_back(ins_time); }
	else {
		while (vsize > 0 && vec[vsize - 1] > ins_time) { vsize--; }
		vec.insert(vec.begin() + vsize, ins_time);
	}
	return;
}

/// FUNCTION TO REMOVE AN EVENT THAT HAS BEEN PASSED
void pregnancy::removetime(vector<double>& rtimevec){
			if(rtimevec.size()>1){
				rtimevec.erase(rtimevec.begin(),rtimevec.begin()+1);
			}
			else{
				rtimevec.clear();
			}
return;
}

///////////////////////////////////////////////// END OF EXPOSURE WITH NO ANC INTERVENTION ////////////////////////


////////////////////////////////////////////////////////// FUNCTIONS FOR ANC INTERVENTIONS //////////////////////
/// FUNCTION TO SEE IF TESTING AND/OR TREATMENT IN FIRST TRIMESTER 
void  pregnancy::any_first_trimester(bool &past) {
	if ((!past)&&((first_tri_rdt > 0) && (gest_time>first_tri_visit_time)&& (first_tri_visit_time >0))) {
		past = true;
		bool current_inf = (start_pcr_time< first_tri_visit_time) && (pcr_clear_time> first_tri_visit_time);
		double first_tri_sens = first_tri_rdt == 1 ? npreg_rdt_sens : perfect_sens;
		if (current_inf&&(runif() < first_tri_sens)) {
			if(ISTeff>runif())drugclearance(first_tri_visit_time, getend_weibull(first_tri_visit_time, ISTscale, ISTshape));
		}
	}
	return;
}

/// FUNCTION TO IMPLEMENT INTERVENTION IN ANC VISITS FROM SECOND TRIMESTER ONWARDS
void pregnancy::IPTISTupdate(int ANC,vector<bool>& past) {
	if ((!past[ANC]) && (gest_time > ANC_times[ANC])) {
		if ((strategy == 0) && (ANC == 1) && ((totalpara > 0) || ((start_pcr_time < ANC_times[ANC]) && (pcr_clear_time > ANC_times[ANC])))) clear_fail = 1;
		if (strategy == 1) {
			// IPTp using "IPT_drug"
			if (IPTeff > runif()) {
				drugclearance(ANC_times[ANC], min(getend_weibull(ANC_times[ANC], IPTscale, IPTshape), gestation_duration - inflag));
			}
			else if(ANC == 1 ){ clear_fail = 1; }
			past[ANC] = true;
		}
		else if (strategy > 1) {
			// INVOLVES AN RDT GET SENSITIVITY
			bool current_infection = ((totalpara > 0 )| ((start_pcr_time< ANC_times[ANC]) & (pcr_clear_time>ANC_times[ANC])));
			double current_sensitivity = perfect_test ? perfect_sens : (ANC != 0 ? (prev_inf ? late_sensitivity_prev_inf[parity] : late_sensitivity[parity]) : sensitivity[histinf]);
			if (strategy == 2) {
				// ISTP - ONLY TREAT POSITIVE WITH "IST_DRUG"
				if (current_infection & (runif() < current_sensitivity)) {
					//DETECTED
					if (ISTeff > runif()) {
						//SUCCESSFULLY TREATED
						drugclearance(ANC_times[ANC], min(getend_weibull(ANC_times[ANC], ISTscale, ISTshape), gestation_duration - inflag));
					}
					else if (ANC == 1) { clear_fail = 1; }
				}
				else if (ANC == 1) { clear_fail = 1; }
			}
			else if (strategy > 2) {
				// HYBRID - STRATEGY 3 IS TEST IN FIRST VISIT IN 2ND TRIMESTER, 4 IS TESTS IN ALL VISITS, ALL WOMEN RECEIVE AT LEAST IPT_DRUG
				if (current_infection & (ANC == 0 || strategy == 4) & (runif() < current_sensitivity)) {
					//DETECTED
					if (ISTeff > runif()) {
						//SUCCESSFULLY TREATED
						drugclearance(ANC_times[ANC], min(getend_weibull(ANC_times[ANC], ISTscale, ISTshape), gestation_duration - inflag));
					}
					else if (ANC == 1) { clear_fail = 1; }
				}
				else {
					if (IPTeff > runif()) {
						//SUCCESSFULLY TREATED
						drugclearance(ANC_times[ANC], min(getend_weibull(ANC_times[ANC], IPTscale, IPTshape), gestation_duration - inflag));
					}
					else if (ANC == 1) { clear_fail = 1; }
				}

			}
		}
	}
	return;
}


/// DRUG CLEARANCE BEGINNING AT CLEARTIME AND PROPHYLAXIS ENDING AT PROPHTIME
void pregnancy::drugclearance(double cleartime, double prophtime) {
	if (gest_time > cleartime) {
		if (start_pcr_time < prophtime) {
			if (start_pcr_time < cleartime) {
				double end = min(cleartime, pcr_clear_time);
				fill(weekperiinf.begin() + (int)start_pcr_time, weekperiinf.begin() + (int)end, 1);
				ever_peri = 1;
			}

			while (start_peri_inftimes_B[0] < prophtime) {
				start_peri_inftimes_B.erase(start_peri_inftimes_B.begin());
				end_peri_time_K.erase(end_peri_time_K.begin());

			}
			if (cleartime < blood_to_IVS) {
				while (plac_exposure_times_P[0] == blood_to_IVS + inflag)plac_exposure_times_P.erase(plac_exposure_times_P.begin());
			}
			else {
				while (plac_exposure_times_P[0] <= prophtime + inflag)plac_exposure_times_P.erase(plac_exposure_times_P.begin());
			}
			pcr_clear_time = 281;
			start_pcr_time = 281;
			gest_time = start_peri_inftimes_B[0];
			next_peri_time = start_peri_inftimes_B[0];
		}

		if ((totalpara > 0) & (startinfdate < prophtime)) {
			end_placenta_active_R.clear();
			start_chronic_placenta.clear();
			totalpara = 0;
			totchron = 0;
			if (startinfdate < cleartime) {
				timeinfected += cleartime - startinfdate;
				fill(weekinf.begin() + (int)startinfdate, weekinf.begin() + (int)cleartime, 1);
				fill(weekperiinf.begin() + (int)startinfdate, weekperiinf.begin() + (int)cleartime, 1);
				ever_peri = 1;
			}
		}
	}
	return;
}


/// GETS WEIBULL DISTRIBUTED DURATION OF PROPHYLAXIS
double pregnancy::getend_weibull(double IPTtime,double scale,double shape){
	return IPTtime+gen_weibull(scale,shape);
}

/// RESETS EVERYTHING READY FOR NEXT PREGNANCY
void pregnancy::clearall(void){
	start_pcr_time = 281;
	pcr_clear_time = 281;
	gest_time=0;
	totalpara=0;
	totchron=0;
	startinfdate=0;
	startchrondate=0;
	next_peri_time=0;
	start_chronic_placenta.clear();
	end_placenta_active_R.clear();
	start_peri_inftimes_B.clear();
	plac_exposure_times_P.clear();
	end_peri_time_K.clear();
	
	proph_fail=0;
	proph_fail_tot=0;
	inf_anc1 = false;
	ever_peri=0;
	ever_plac=0;
	timechron=0;
	timeinfected=0;
	clear_fail = 0;
	total_plac = 0;

	fail = 0;
	prev_inf=0;
	fill(weekinf.begin(),weekinf.end(),0);
	fill(weekperiinf.begin(),weekperiinf.end(),0);
return;
}

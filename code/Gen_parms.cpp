


void gen_parms::setup(void){
	/// INITIALISE PARAMETERS ORDERED AS IN TABLE S1
	/// AGE AND HETEROGENEITY
	numhet = 5;
	agebracks = 36;
	// set up heterogeneity approximation
	rho = 0.85;
	a0 = 8.0 * 365.0;
	sigma2 = 1.03885;
	//INFECTIOUS PERIODS
	rd = 1.00 / 5.00;
	ra = 1.00 / 195;
	ru = 1.00 / 107.867;
	rt = 1.00 / 5.00;
	rp = 1.00/9.00;

	//IMMUNITY PARAMETERS
	beta=0.930766;
	bmin=0.00640231;
	Ib0 = 1072.1427;
	tau_b = 0.20554003;
	wb = 3650;
	Kb = 6.0487;
	d_min = 160;
	
	IA0 = 4732.5;
	Ka = 5.0;
	wa=3650;
	
	Ic0 = 53.03;
	Kc = 2.01471;
	wca = 10950;
	pcm = 0.532052;
	wcm = 230.402;
	
	///EXPONENTIAL DEMOGRAPHY WITH 21 YEAR MEAN AGE
	DY = 365;
	mu=1/(21*DY);	
	
	
	
	//RESIZE STATE SPACE VECTORS
	FOI.resize(agebracks, vector<double>(numhet));
	den_het = S = A = U = D = P = T = phis = IB = FOI;
	xI.resize(agebracks);
	recI.resize(agebracks);
	het_x.resize(numhet);
	het_wt.resize(numhet);
	
	/// DEFINE GAUSSIAN HERMITE WEIGHTS
	het_x[0] = -2.85697001387;
	het_x[1] = -1.35562617997;
	het_x[2] = 0;
	het_x[3] = -het_x[1];
	het_x[4] = -het_x[0];
	het_wt[0] = 0.0112574113277;
	het_wt[1] = 0.222075922006;
	het_wt[2] = 0.533333333333;
	het_wt[3] = het_wt[1];
	het_wt[4] = het_wt[0];
	
	///DEMOGRAPHY AND AGE-DEPENDENT EXPOSURE SET-UP
	deltas.push_back(mu);
	double sumage=0;
	double pop=0;
		for(int i=0; i<agebracks;i++){
		ageFOI.push_back(1-rho*exp(-sumage/a0));
				if(i<8){
		age_width.push_back(365*0.25);
		}
		else if(i<16){
		age_width.push_back(365*0.5);
		}
		else if(i<25){
		age_width.push_back(365);
		}
		else if(i<agebracks){
		age_width.push_back(5*365);
				}
		sumage+=age_width[i];
		if(i==0){
			demog.push_back(1/(1+1/(age_width[i]*mu)));
				}
		if(i<agebracks-1){
		age_rate.push_back(1/age_width[i]);
			if(i>0){
		deltas.push_back(age_rate[i-1]);
		demog.push_back(demog[i-1]*age_rate[i-1]/(age_rate[i]+mu));
		}
		gammas.push_back(mu+age_rate[i]);
	}
	else{
		age_rate.push_back(0);
	gammas.push_back(mu);
	deltas.push_back(age_rate[i-1]);
	demog.push_back(demog[i-1]*age_rate[i-1]/(age_rate[i]+mu));
	}	
	for(int j=0;j<numhet;j++){
			den_het[i][j]=demog[i]*het_wt[j];
			pop+=den_het[i][j];
	}
	
	}
		return;
 }

void gen_parms::init_determ_model(void){
/// INITIALISE PARAMETERS ORDERED AS IN TABLE S1
	vector<double> rel_foi(numhet);
	for(int i=0;i<numhet;i++){
		rel_foi[i]=exp(-sigma2/2+pow(sigma2,0.5)*het_x[i]);
	}
	double EIRy=EIR/365;
///// GET EQUILIBRIUM EIR AND IMMUNITY FUNCTIONS
	vector<vector<double>> eqEIR(agebracks, vector<double>(numhet));
	vector<vector<double>> Ib(agebracks, vector<double>(numhet));
	vector<vector<double>> Ic(agebracks, vector<double>(numhet));
	vector<vector<double>> Ica(agebracks, vector<double>(numhet));
	vector<vector<double>> Icm(agebracks, vector<double>(numhet));
	vector<vector<double>> nondis(agebracks, vector<double>(numhet));
	vector<double> Ia(agebracks);
	double rmaxmin=1/(ra*d_min);
	double denagetrack=1.0;
	double sumage=0.0;
	for(int i=0; i<agebracks;i++){
		Ia[i]= wa *(1-exp(-1.0*sumage/wa));
		double Bli;
		if(Ia[i]==0){
			Bli =0;
		}
		else{
			Bli =pow(Ia[i]/IA0,Ka);
		}
		recI[i]=ra*(1+(rmaxmin -1)*Bli/(1+Bli));
		//cout << recI[i] <<"  "<< Bli<<"\n";
		sumage+=1.0*age_width[i];
		if(i==0){
		xI[i]=demog[i]/mu;
		}
		else{
			denagetrack=demog[i-1]*age_rate[i-1];
			xI[i]=demog[i]/denagetrack;
		}
		for(int j=0;j<numhet;j++){
			eqEIR[i][j]=EIRy*ageFOI[i]*rel_foi[j];
			if(i==0){
				Ib[i][j]=pow(eqEIR[i][j],tau_b)*xI[i]/(1+xI[i]/wb);
			}
			else{
				Ib[i][j]=(Ib[i-1][j]+pow(eqEIR[i][j],tau_b)*xI[i])/(1+xI[i]/wb);
			}
			FOI[i][j] = eqEIR[i][j] * ((beta - bmin) / (1 + pow(Ib[i][j] / Ib0, Kb)) + bmin);
			if (i == 0) {
				Ica[i][j] = FOI[i][j] * xI[i] / (1 + xI[i] / wca);
			}
			else {
				Ica[i][j] = (Ica[i - 1][j] + FOI[i][j] * xI[i]) / (1 + xI[i] / wca);
			}
		}
	}
	for(int i=0; i<agebracks;i++){
		for(int j=0;j<numhet;j++){
			if(i==0){
				Icm[i][j]=pcm* Ica[27][j]/(1+xI[i]/wcm);
			}
			else{
				Icm[i][j]= Icm[i-1][j]/(1+xI[i]/ wcm);
			}
			Ic[i][j]= Icm[i][j]+ Ica[i][j];
			phis[i][j]=1/(1+pow(Ic[i][j]/Ic0,Kc));
		}
	}
	///// GET EQUILIBRIUM STATES
	vector<vector<double>> propd(agebracks,vector<double> (numhet));
	vector<vector<double>> propt(agebracks,vector<double> (numhet));
	vector<vector<double>> propp(agebracks,vector<double> (numhet));
	double pop=0;
	for(int i=0;i<agebracks;i++){
		for(int j=0;j<numhet;j++){
			propt[i][j]=FOI[i][j]*phis[i][j]*(ft)/(gammas[i]+rd);
			propd[i][j]= FOI[i][j]*phis[i][j]*(1-ft)/(gammas[i]+rd);
			propp[i][j]=rd*propt[i][j]/(gammas[i]+rp);
			if(i==0){
				nondis[i][j]=den_het[i][j]/(1+propd[i][j]+propt[i][j]+propp[i][j]);
				D[i][j]=propd[i][j]*nondis[i][j];
				T[i][j]=propt[i][j]*nondis[i][j];
				P[i][j]=propp[i][j]*nondis[i][j];
				A[i][j]=(nondis[i][j]*(1-phis[i][j])*FOI[i][j] +rd*D[i][j])/(FOI[i][j]+recI[i]+gammas[i]);
				U[i][j]=recI[i]*A[i][j]/(FOI[i][j]+ru+gammas[i]);
				S[i][j]=nondis[i][j]-A[i][j]-U[i][j];
			}
			else{
				nondis[i][j]=(den_het[i][j]-deltas[i]*(D[i-1][j]/(rd+gammas[i])+T[i-1][j]/(rd+gammas[i])+(rd*T[i-1][j]/(rd+gammas[i])+P[i-1][j])/(rp+gammas[i])))/(1+propd[i][j]+propt[i][j]+propp[i][j]);
				T[i][j]=propt[i][j]*nondis[i][j]+deltas[i]*(T[i-1][j]/(gammas[i]+rd));
				D[i][j]=propd[i][j]*nondis[i][j]+deltas[i]*(D[i-1][j]/(gammas[i]+rd));	
				P[i][j]=propp[i][j]*nondis[i][j]+deltas[i]*(rd*T[i-1][j]/(gammas[i]+rd)+P[i-1][j])/(gammas[i]+rp);	
				A[i][j]=(deltas[i]*A[i-1][j]+nondis[i][j]*(1-phis[i][j])*FOI[i][j]+rd*D[i][j])/(FOI[i][j]+recI[i]+gammas[i]);
				U[i][j]=(deltas[i]*U[i-1][j]+recI[i]*A[i][j])/(FOI[i][j]+ru+gammas[i]);
				S[i][j]=nondis[i][j]-A[i][j]-U[i][j];
			}
		}
	}
	///// CALCULATE PCR PREVALENCES
	double prop_adult=0;
	pcr_adult=0;
	for(int i=25;i<32;i++){
	for(int j=0;j<numhet;j++){
		pcr_adult+=(A[i][j]+U[i][j]+D[i][j]+T[i][j]);
		}
		prop_adult+=demog[i];
	}
	pcr_adult/=prop_adult;
	double child=0;
	pcr_child=0;
	for(int i=0;i<12;i++){
		for(int j=0;j<numhet;j++){
			pcr_child+=(A[i][j]+U[i][j]+D[i][j]+T[i][j]);
		}
		child+=demog[i];
	}
	pcr_child/=child;
	return;
}


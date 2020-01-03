//// STRUCTURE DETERMINING GENERAL POPULATION MODEL ////


struct gen_parms{
	
	double rho;
	double a0;
	double sigma2;
	
	double rd;
	double ra;
	double ru;
	double rt;
	double rp;

	//IMMUNITY PARAMETERS
	double beta ;
	double bmin;
	double Ib0;
	double tau_b;
	double wb;
	double Kb;
	double d_min;

	double IA0;
	double Ka;
	double wa;

	double Ic0;
	double Kc;
	double wca;
	double pcm;
	double wcm;

	double mu;
	

 double DY;
 
 int numhet;
 vector<double> het_x;
 vector<double>  het_wt;
 vector<vector<double>>  den_het;

 int agebracks;
 vector<double> age_width;
 vector<double> age_rate;
 vector<double> gammas;
 vector<double> deltas;
 vector<double> demog;

 //int sites;
 double EIR;
 double ft;


 double pcr_adult;
 double pcr_child;

  ///----Initial Population Prevs ---/////
	vector<vector <double>> S,A,U,D,P,T,phis,IB,FOI;
	vector <double> xI,recI,ageFOI;	

  	 
 void setup(void);
 
	
 void init_determ_model(void);
 

};



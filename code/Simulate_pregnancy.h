
/// FOR RUNNING PREGNANCIES ///

struct pregnancy {
	double gest_time;
	double gestation_duration;
	double blood_to_IVS;
	/// RELEVANT NON-PREG SPECIFIC PARAMS
	double age;
	int age_cat;
	double parity;
	int hetcat;
	int start_inf_state;
	double pcr_prev;
	double cFOI;
	double recA;
	double rd;
	double ru;
	double rp;
	double ft;
	double deathprob;
	///-----------------PREG-SPECIFIC IMMUNITY PARAMETERS-----------------------////
	int histplac;
	int histinf;
	double inflag;
	double prog1;
	double prog2;
	double Kp;
	double Ip0;
	double Kchron;
	double Ichron0;
	double fracD;
	double rate;
	
	///-----------------NATURAL HISTORY OF PREGNANCY-----------------------////

	double fracA[200];
	double progC[200];

	double p1;
	double p2;

	double totalpara;
	double totchron;




	double ever_peri;
	double ever_plac;
	double total_plac;
	double timechron;
	double timeinfected;
	double timeperi;
	double clear_fail;
	
	double proph_fail;
	double proph_fail_tot;
	vector<double> weekinf;
	vector<double> weekperiinf;

	///// PERIPHERAL INFECTION
	//start and end of peripheral infections (B and K in supplementary info)
	vector<double> start_peri_inftimes_B;
	vector<double> end_peri_time_K;
	//for tracking peripheral events
	double next_peri_time;
	double start_pcr_time;
	double pcr_clear_time;

	/// PLACENTAL INFECTION
	//start and end of placental infections (P and R in supplementary info)
	int init_plac;
	vector<double> plac_exposure_times_P;
	vector<double> end_placenta_active_R;
	// beginning of chronic infection
	vector<double> start_chronic_placenta;
	//for tracking placental infection
	double nextplactime;
	int nextplacchange;
	double startinfdate;
	double startchrondate;

	
	double fail;
	double failprob;
	bool prev_inf;

///// NATURAL PROGRESSION FUNCTIONS  /////
void generateperiinfs(void);
void run_to_delivery(void);
double nextstate(void);
void removetime(vector<double>& rtimevec);
void newinfpath(void);
void placentalevent(void);


void placevec(vector<double> &vec, double ins_time);
void clearall(void);

///// INTERVENTION PARAMETERS  /////
/* strategy to implement
0: No intervention
1: IPTp (with IPTp drug)
2: ISTp (with ISTp drug)
3: Hybrid SSTp (at first ANC test and treat with an ISTp drug if test-positive, IPTp drug otherwise and at all other visits)
4: Hybrid ISTp (test at all visits, ISTp drug if test-positive, IPTp drug otherwise)
*/
int strategy;
/// Gestation time of each ANC visit
vector<double> ANC_times;
// whether intervention disrupts immunity (or treating intervention as newly implemented across gravidities
bool intervene_imm;
// a counter if immunity is not dirupted (or treating intervention as newly-implemented
int intervene;
/// whether to implement first trimester testing if an early visit
int first_tri_rdt;
/// time of first trimester visit
int first_tri_visit;
double first_tri_visit_time;

double npreg_rdt_sens;
int perfect_test;
double perfect_sens;
double IPTeff;
double IPTscale;
double IPTshape;
double ISTeff;
double ISTscale;
double ISTshape;


void any_first_trimester(bool &past);
void IPTISTupdate(int ANC,vector<bool> & past);

double primi_sens;
double rdt_preg_offset;
double rdt_preg_shape;
double rdt_adult_offset;
double rdt_adult_shape;
double late_PG_sens_odds;
double late_grav_OR;
double prev_inf_OR;
vector<double> sensitivity;
vector<double> late_sensitivity;
vector<double> late_sensitivity_prev_inf;
void drugclearance(double IPTtime,double prophtime);
double getend_weibull(double IPTtime, double scale, double shape);

//// FOR OUTPUTS



};



#include "standard.h"
#ifndef STANDARD_C
#define STANDARD_C

string extract_time_info(const char* format){
	char c[40];
	time_t tim=time(0);
	strftime(c, 40, format, localtime(&tim));
	return string(c);
}
string timestr(void){
	return extract_time_info("%H.%M.%S");
}
string datestr(void){
	return extract_time_info("%A %d %b %Y");
}


void errormsg(const string s){
	cout<<s<<endl;
	system("PAUSE");
	exit(1);
	}
void msg(const string s){
	cout<<s<<endl;
	}
template <class T>
string getname(const T &name){
	stringstream s;
	s<<name;
	return s.str();}
string as_string(const size_t &k){
	stringstream s;
	s << static_cast<unsigned int>(k);
	return s.str();
}
void error_crit(const string s){
	cout << '\n' << s << endl;
	system("PAUSE");
	exit(1);
}
void error_crit(void){
	error_crit("Critical error");
}

///- Rounds double to nearest integer --///
int round2(double z){
	return int (z+0.5);}




///--- Beginning of Random number generators taken from Numerical Recipes Book (runif=uniform(0,1) etc)---/////
long SEED;
double powgam;
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 3.0e-16
#define RNMX (1.0-EPS)
double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


double runif() // Generates a uniform rv ~U(0,1)
{return ran1(&SEED);}


double bingen(double n, long double p){// Generates a binomial rv ~Bin(n,p)
	double x=0;
	for(double k=0;k<n;k++){
		if(runif()<p){
			x++;}
	}
	return x;
}
double genexp(double lambda){
	return -log(runif())/lambda;
}
double gengam(int shape, double scale){
	double x=0;
	for(int i=0;i<shape;i++){
		x+=-log(runif());
	}
	return scale*x;
}
double gen_weibull(double scale,double shape){
	return(pow(-1/scale*(log(1-runif())),1.0/shape));}

///--- End of Random number generators ---/////



#endif



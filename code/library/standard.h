#ifndef STANDARD_H
#define STANDARD_H

#include <random>
#include <math.h> 
#include <cmath>
#include <vector>
#include <queue>
#include <deque>
#include <list>
#include <cstdio>
#include <cstdlib>
#include <ios>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <utility>
#include <assert.h>
#include <iterator>
#include <numeric>
#include<sstream>
#include <strstream>


using namespace std;
void errormsg(const string s);
void msg(const string s);
template <class T>
string getname(const T &name);
string as_string(const size_t &k);
template <class T>
string as_string(const T &k){
	stringstream s;
	s << k;
	return s.str();
}
void error_crit(const string s);
void error_crit(void);

int round2(double z);

double runif();
double bingen(double n, long double p);
double genexp(double lambda);
double gengam(int shape, double scale);
double gen_weibull(double scale,double shape);
string extract_time_info(const char* format);
string timestr(void);
string datestr(void);

#endif
#ifndef _MY_GLOBAL_
#define _MY_GLOBAL_

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>

//#define FERMION
#define SYMMETRY 3

#define TOL 1e-10
#define NEV 6
#define ZERO 0.01  // for print only

using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;
using std::vector;
using std::complex;
using std::string;
using std::abs;

extern string filename;
extern bool print;
extern int symmetry_sector;

#endif

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

//#define PBC
#define FERMION
#define SYMMETRY 2
//#define S4E

#define TOL 1e-10  // the tolerence for arpack solver
#define ACC 10  // the precision of data printed
#define NEV 6  // the number of eig-vals to be claculated/printed
#define NSI 20  // the number of singular values to be printed
#define ZERO 0.01  // for print only
#define PRINT_TO_SCREEN 0
#define PI 3.14159265359

using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;
using std::ifstream;
using std::fstream;
using std::ios;
using std::vector;
using std::complex;
using std::string;
using std::abs;

#endif

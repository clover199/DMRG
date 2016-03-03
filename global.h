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

#define PBC
#define FERMION
#define SYMMETRY 2

#define TOL 1e-10
#define ACC 10  // the precision of data printed
#define NEV 6  // the number of eig-vals to be claculated/printed
#define NSI 10  // the number of singular values to be printed
#define ZERO 0.01  // for print only
#define PI 3.14159265359

using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;
using std::vector;
using std::complex;
using std::string;
using std::abs;

#endif

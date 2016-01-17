#ifndef _MY_FUNCTIONS_
#define _MY_FUNCTIONS_

#include "global.h"
#include "mps.h"
#include "mpo.h"

// We assume there are at most 8 tunable parameters in the Hamiltonian.
// The function puts the names of the parameters into a vector.
vector<string> set_para_name(string s1="", string s2="", string s3="", string s4="",
                             string s5="", string s6="", string s7="", string s8="");


// The function gets the parameters as the same order of the parameter names
vector<double> set_para_val(int argc_, char *argv_[], int& sites, int& cutoff, int& sweep,
                            int& sector, vector<string>& para_name);


void update_two(int l, int r, int cutoff, mps<double>& my_mps, mpo<double>& my_mpo);


void dmrg(mps<double> my_mps, mpo<double> my_mpo, int cutoff, int sweep);

#endif

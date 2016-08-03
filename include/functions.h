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
                            int& sector, vector<string>& para_name, string& filename);

template <typename T>
qtensor<T> combine_mpo(qtensor<T>& H1, qtensor<T>& H2);

template <typename T>
void ed(mps<T>& my_mps, mpo<T>& my_mpo,
        int cut, int swp, int sector, const string& filename);
/*
template <typename T>
void init_dmrg1(mps<T>& my_mps, mpo<T>& my_mpo,
                int cutoff, int sweep, int sector, const string& filename);

template <typename T>
void init_dmrg2(mps<T>& my_mps, mpo<T>& my_mpo,
                int cutoff, int sweep, int sector, const string& filename);*/

template <typename T>
void dmrg(mps<T>& my_mps, mpo<T>& my_mpo,
          int cutoff, int sweep, int sector, const string& filename);

template <typename T>
void dmrg2(mps<T>& my_mps, mpo<T>& my_mpo,
           int cutoff, int sweep, int sector, const string& filename);

template <typename T>
void idmrg(mps<T>& my_mps, mpo<T>& my_mpo,
          int cutoff, int sweep, int sector, const string& filename);

template <typename T>
void idmrg1(mps<T>& my_mps, mpo<T>& my_mpo,
           int cutoff, int sweep, int sector, const string& filename);

template <typename T>
T calc(mps<T>& my_mps);

template <typename T>
T calc(mps<T>& l_mps, mps<T>& r_mpo);

template <typename T>
T calc(mps<T>& my_mps, mpo<T>& my_mpo);

#endif

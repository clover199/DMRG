#ifndef _MY_USEFUL_
#define _MY_USEFUL_

#include <vector>
using std::vector;

#include "qtensor.h"

// define the rule for adding symmetries
int ss(int a, int b);

/* define the fermion parity of each symmetry sector */ 
int ff(int a);

void print_matrix(const vector< vector<int> >& map);

void print_energy(ofstream& data, int l ,int r, double* val, int n=NEV);

void print_singular(ofstream& data, int l, int r, const vector<double>& s, int cutoff);

void generate_map(vector< vector<int> >& map, const vector<int>& index);

int get_dimension(const vector< vector<int> >& dim,
                  const vector< vector<int> >& sym);

template <typename T>
void generate_dim_sym(vector< vector<int> >& dim, vector< vector<int> >& sym,
                      const qtensor<T>& lenv, int lnum,
                      const qtensor<T>& renv, int rnum, int sector);

#endif

#ifndef _MY_QTENSOR_
#define _MY_QTENSOR_

#include "global.h"
#include "tensor.h"

// Attention: the constructor function and functions:
// 'exchange', 'transpose', 'cut'
// can only accept at most 8 indexes

// T = 'double', 'complex<double>'

template <typename T> 
class qtensor
{
  template <typename> friend class qtensor;

private:
  vector< vector<int> > dim_;  // the dimension for different sectors
                               // dim_[i][j] is the dimension for index i, symmetry j
  vector< vector<int> > sym_;  // symmetries are represented by 0,1,2...
  vector< tensor<T> > val_;
  
  // get/check "dim_" from "sym_" and "val_"
  void check_dim_(int n);
#ifdef FERMION
  vector<int> dir_; // -1 means left, 1 means right, 0 meas not a physical index
#endif

public:
  // claim space with all elements as zero
  qtensor(int i0=0, int i1=0, int i2=0, int i3=0, 
          int i4=0, int i5=0, int i6=0, int i7=0);

  qtensor(T* val, const vector< vector<int> >& dim,
          const vector< vector<int> >& sym);

  ~qtensor() {
    dim_.clear();
    sym_.clear();
    val_.clear();
#ifdef FERMION
    dir_.clear();
#endif
  }

  void update(tensor<T>& val,
              int i0, int i1=-1, int i2=-1, int i3=-1, 
              int i4=-1, int i5=-1, int i6=-1, int i7=-1);

  void clear() {
    dim_.clear();
    sym_.clear();
    val_.clear();
  }

  // only used for the 4-tensor MPO
  qtensor id() const;

  // only used for the 3-tensor MPO
  qtensor id(const qtensor<T>& bulk_mpo) const;

#ifdef FERMION
  void add_sign(int i0=2, int i1=2, int i2=2, int i3=2, 
                int i4=2, int i5=2, int i6=2, int i7=2);
#endif
  
  vector< vector<int> > dimension_all() const {return dim_;}
  
  int dimension(int n) const {
    int ret = 0;
    for(int i=0;i<dim_[n].size();i++) ret += dim_[n][i];
    return ret;
  }
  
  int dimension() const {
    int ret = 1;
    for(int i=0;i<dim_.size();i++) ret *= dimension(i);
    return ret;
  }

  int val(T* out, vector< vector<int> >& dim, vector< vector<int> >& sym) const {
    int k = 0;
    for(int i=0;i<val_.size();i++) for(int j=0;j<val_[i].val_.size();j++)
       out[k++] = val_[i].val_[j];
    dim = dim_;
    sym = sym_;
    return k;
  }
  
  bool empty() const {
    if(val_.size()==0) return true;
    else return false;
  }
  
  void print() const;
  
  // all=false: print only non-zero values
  // all=true:  print all values including zeros
  void print_all(bool all=false) const;

  void print_matrix() const;

  qtensor conjugate();

  // conver the double qtensor to complex
  qtensor<complex<double> > comp() const;  

  // remove the all-zero symmetry sectors
  qtensor simplify() const;

  // exchange the two indexes at positions a, b (start from 0)
  // for a two-tensor and a=0, b=1, this is just matrix transpose.
  qtensor exchange(int a=0, int b=1) const;

  // transpose as a matrix with the first num indexes as left index
  // and the rest indexes as right index
  qtensor shift(int num=1) const;

  // take the kth elements of index n
  qtensor take(int n, int k) const;

  // only used for the 4-tensor MPO
  qtensor left() const;

  // only used for the 4-tensor MPO
  qtensor right() const;
    
  // combine the index from min to max (including max)
  // can only combine indexes with the same (or no) symmetries 
  qtensor<T> combine(int min=0, int max=1) const;
  
  // index is the index to be splitted
  // dim is the dimensions for the splitted part
  // sym is the sym for all the index
  qtensor<T> split(int index, vector< vector<int> > dim) const;

  qtensor<T> remove_symmetry() const;

  // retain 'cut' number of index 'index'
  qtensor cut(int index, vector<int> cutoff) const;
  
  // result = alpha * A + beta * B
  qtensor& plus(const qtensor& A, const qtensor& B, T alpha=1, T beta=1);

  T trace(qtensor& A, bool fermion=false);
  
  qtensor operator+(const qtensor& A) const {
    qtensor<T> ret;
    ret.plus(*this, A);
    return ret;
  }
  
  // contract tensor A and B with the left/right most index
  // This function cannot be used for fermion operator contraction
  qtensor& contract(qtensor& A, qtensor& B,
                    char transa='N', char transb='N', int num=1);
  
  qtensor operator*(qtensor& A) {
    qtensor<T> ret;
    ret.contract(*this, A);
    return ret;
  }
  
  // contract index a of tensor A with index b of tensor B
  // change is only used for fermions. change=true: exchange A,B position
  qtensor& contract(const qtensor& A, int a, const qtensor& B, int b,
                    bool change = false);
  
  // multiply the tensor by an array (as a matrix)
  // used in the sparse matrix multiplication routine in dmrg 
  // map can be generated from function get_map
  // map: val_ index, in begin, out begin, row, col, beta
  void contract(T* out, T* in,
                const vector< vector<int> >& map,
                char transa='N', char transb='N', int num=1);

  bool get_map(vector< vector<int> >& map_ret, 
               vector< vector<int> >& dim_ret, vector< vector<int> >& sym_ret,
               const vector< vector<int> >& dim, const vector< vector<int> >& sym,
               char transa='N', char transb='N',
               int num=1, double beta=0, bool back=false) const;

  // SVD between first n and the rest indexes
  // A = U * S * V
  // returns the number of singular values (length of S)
  vector<double> svd(qtensor& U, qtensor<double>& S, qtensor& V,
                     int num=1, int cutoff=0);

  // calculate the eigenvalues and the eigenvector for the lowest eig-val
  int eig(double* val, T* vec, int sector=-1);
};
#endif

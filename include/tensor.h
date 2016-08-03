#ifndef _MY_TENSOR_
#define _MY_TENSOR_

#include "global.h"

/* Attention: the constructor function and functions:
    'update', 'id'
   can only accept at most 8 indexes 

   A q-tensor has the form:
       A[i0,i1,i2,...,iq]
   with i0, i1, i2, ..., iq as its index number.
   
   T = 'double', 'complex<double>' */

template <typename T> 
class tensor
{
  template <typename> friend class tensor;
  template <typename> friend class qtensor;

private:

  vector<int> index_;  /* dimension of each index */
  vector< T > val_;  /* values of the tensor rearranged as an array*/
  
  inline T* begin() { return &*val_.begin(); }
  
  int dim_() const {
    if(index_.size()==0) return 0;
    int ret = 1;
    for(int i=0;i<index_.size();i++) ret *= index_[i];
    return ret;
  }

  int loc_(const vector<int>& in) const {
    if(in.size()!=index_.size())
    {
      cerr << "Error in tensor loc_: wrong dimension." << endl;
      return 0;
    }
    int ret = 0;
    for(int i=0;i<in.size();i++) ret = ret*index_[i]+in[i];
    return ret;
  }

  void refine_(vector<int>& loc) const {
    if(loc.size()!=index_.size())
    {
      cerr << "Error in tensor loc_: wrong dimension." << endl;
      return ;
    }
    for(int i=index_.size()-1;i>0;i--) if(loc[i]>=index_[i])
    {
      loc[i-1] += loc[i] / index_[i];
      loc[i] = loc[i] % index_[i];
    }
    if(loc.size()) loc[0] = loc[0] % index_[0];
  }

public:

  /* claim space with all elements as zero */
  tensor(int i0=0, int i1=0, int i2=0, int i3=0, 
         int i4=0, int i5=0, int i6=0, int i7=0);

  tensor(const vector<int>& index);
  
 /* create tensor with input values */
  tensor(T* val,
         int i0,   int i1=0, int i2=0, int i3=0, 
         int i4=0, int i5=0, int i6=0, int i7=0);
		 
  tensor(T* val, const vector<int>& index);

  ~tensor() { index_.clear(); val_.clear(); }

  /* assign value to a particular member */ 
  void update(T val,
        int i0,    int i1=-1, int i2=-1, int i3=-1, 
        int i4=-1, int i5=-1, int i6=-1, int i7=-1);

  void update(T val, const vector<int>& loc);

  /* assign zero/one/random values to the tensor */
  void zeros();
  void ones();
  void rand();

  /* set A=1 only for index satisfying a1=b1, a2=b2, ...
     if no arguments given, set A=1 when all indexes are the same */
  void id(int a1=-1, int b1=-1, int a2=-1, int b2=-1,
          int a3=-1, int b3=-1, int a4=-1, int b4=-1);

/* *********** print and get properties of the tensors ********** */

  int index() const {return index_.size();}

  int dimension(int n) const {return index_.at(n);}
  
  int dimension() const {return val_.size();}

  T val(int i) const {return val_.at(i);}

  T val(const vector<int> loc) const {return val_.at( loc_(loc) );}
  
  /* test if all the values are zero */
  bool empty() const;

  void print() const;
  
  void print_matrix() const;

  void print_index() const; 

  void write_file(fstream& name) const;

  void read_file(fstream& name);

  void write_file(ofstream& name) const;

  void read_file(ifstream& name);

/* *********** calculation and changes to the tensors ********** */

  void replace(T val, int i) { val_.at(i) = val; }

  void replace(T val, const vector<int>& loc) { val_.at( loc_(loc) ) = val; }

  void conjugate();

  void to_double(tensor<double>& ret) const;

  void to_complex(tensor<complex<double> >& ret) const;

  /* exchange the two indexes at positions a, b (start from 0) for
     a two-tensor and a=0, b=1, this is just matrix transpose. */
  void exchange(tensor& ret, int a=0, int b=1) const;

  /* transpose as a matrix with the first num indexes as left
     indexes and the rest indexes as right index */
  void shift(tensor& ret, int num=1) const;
  
  /* take the k th elements of index n */
  void take(tensor& ret, int n, int k) const;
  
  /* turn to a num-tensor with only diagonal elements non-zero 
     the tensor has dimension min for all indexes */
  void diag(int num=1);
  
  /* combine the index from min to max (including min and max) */
  void combine(int min=0, int max=1);
  
  /* resize index k to be size n
     cut the values (small n) or add zeros (large n) */
  void resize(int k, int n);
  
  /* resize all the index and set all values to be 0 */
  void resize(const vector<int>& index);

  /* returns sum_i val[i,i,...,i] */
  T trace() const;
  
  /* returns sum val[i]*A[i] for all elements */
  T trace(tensor& A) const;
  
  /* returns: alpha * A + beta * B */
  void plus(tensor& ret, const tensor& B, T alpha=1., T beta=1.) const;
  
  tensor operator+(const tensor& A) const {
    tensor<T> ret;
    plus(ret, A);
    return ret;
  }
  
  tensor operator-(const tensor& A) const {
    tensor<T> ret;
    plus(ret, A, 1.0, -1.0);
    return ret;
  }

  /* returns alpha * A + beta */
  void times(T alpha, T beta=0.0);
  
  tensor operator*(T alpha) const {
    tensor<T> ret = *this;
    ret.times(alpha);
    return ret;
  }
  
  tensor operator/(T alpha) const {
    tensor<T> ret = *this;
    ret.times(1./alpha);
    return ret;
  }
  
  void gemm(T* out, T* in, char transa, char transb, 
            int m, int n, int k, T alpha, T beta);

  /* contract tensor A and B with the num left/right most index
     if cover = false: the result is added to the original tensor
     if cover = true: the result covers the original tensor. */
  void contract(tensor& A, tensor& B, char transa='N', char transb='N',
                   int num=1, bool cover=false);
  
  tensor operator*(tensor& A) {
    tensor<T> ret;
    ret.contract(*this, A, 'N', 'N', 1, true);
    return ret;
  }
  
  /* contract index a of tensor A with index b of tensor B
     returns: A*B*sign
     the final index has the form:
         A0 A1 A2 ... Aa-1 B0 B1 ... Bm Aa+1 ... An
     if cover = false, the result tensor will be added to the original one
     if cover = true, the original tensor will be destroyed. */
  void contract(const tensor& A, int a, const tensor& B, int b,
                   bool cover=false, double sign = 1);
  
  /* multiply the tensor A by an array 'in' (as a matrix)
     only multiply the 'num' left/right most index
     for the array, use the same notation as matrix multiplication
     returns:  out = alpha * A * in + beta * out  */
  void contract(T* out, T* in, int row, int col,
                char transa='N', char transb='N', int num=1,
                double beta=0, double alpha=1);

  void gesvd(T* in, int row, int col, double* s, T* u, T* v);

  /* SVD between first num and the rest indexes
         A = U * S * V
     returns the number of singular values (length of S) */
  int svd(tensor& U, tensor<double>& S, tensor& V, int num=1);

  void heev(int n, double* val, T* vec);  

  int eig(double* val, T* vec);

};

#endif

#ifndef _MY_TENSOR_
#define _MY_TENSOR_

#include "global.h"

/* Attention: the constructor function and functions:
     'exchange', 'transpose', 'cut'
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
    if(in.size()!=index_.size()) return 0;
    int ret = 0;
    for(int i=0;i<in.size();i++) ret = ret*index_[i]+in[i];
    return ret;
  }
public:
  /* claim space with all elements as zero */
  tensor(int i0=0, int i1=0, int i2=0, int i3=0, 
         int i4=0, int i5=0, int i6=0, int i7=0);
		 
  /* create tensor with input values */
  tensor(T* val,
         int i0,   int i1=0, int i2=0, int i3=0, 
         int i4=0, int i5=0, int i6=0, int i7=0);

  ~tensor() {
    index_.clear();
    val_.clear();
  }

  /* assign value to a particular member */ 
  void update(T val,
        int i0,    int i1=-1, int i2=-1, int i3=-1, 
        int i4=-1, int i5=-1, int i6=-1, int i7=-1);

  /* assign random values to the tensor */
  void rand();

  int index() const {return index_.size();}

  int dimension(int n) const {return index_[n];}
  
  int dimension() const {return val_.size();}

  T val(int i) const {return val_[i];}

  void print(bool all=false) const;
  
  void print_matrix() const;

  void print_index() const {
    if(index_.size()==0) cout << "()" << endl;
    else
    {
      cout << "(" << index_[0];
      for(int i=1;i<index_.size();i++) cout << ", " << index_[i];
      cout << ")\n";
    }
  }

  void write_file(fstream& name) const;

  void read_file(fstream& name);

  void write_file(ofstream& name) const;

  void read_file(ifstream& name);

  /* changes the original tensor */
  tensor& conjugate();

  /* convert the double tensor to complex */
  tensor<complex<double> > comp() const;

  /* exchange the two indexes at positions a, b (start from 0) for
     a two-tensor and a=0, b=1, this is just matrix transpose. */
  tensor exchange(int a=0, int b=1) const;

  /* transpose as a matrix with the first num indexes as left
     indexes and the rest indexes as right index */
  tensor shift(int num=1) const;
  
  /* turn to a num-tensor with only diagonal elements non-zero 
     the tensor has dimension min for all indexes */
  tensor& diag(int num=1);
  
  /* combine the index from min to max (including min and max) */
  tensor& combine(int min=0, int max=1);
  
  /* resize index k to be size n */
  tensor resize(int k, int n) const;
  
  /* returns: alpha * A + beta * B */
  tensor& plus(const tensor& A, const tensor& B, T alpha=1, T beta=1);
  
  tensor operator+(const tensor& A) const {
    tensor<T> ret;
    ret.plus(*this, A);
    return ret;
  }

  T trace(tensor& A);

  tensor times(T alpha) const {
    tensor<T> ret = *this;
    for(int i=0;i<val_.size();i++) ret.val_[i] *= alpha;
    return ret;
  }
  
  /* contract tensor A and B with the num left/right most index
     if cover = false: the result is added to the original tensor
     if cover = true: the result covers the original tensor. */
  tensor& contract(tensor& A, tensor& B, char transa='N', char transb='N',
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
  tensor& contract(const tensor& A, int a, const tensor& B, int b,
                   bool cover=false, double sign = 1);
  
  /* multiply the tensor A by an array 'in' (as a matrix)
     only multiply the 'num' left/right most index
     for the array, use the same notation as matrix multiplication
     returns:  out = alpha * A * in + beta * out  */
  void contract(T* out, T* in, int row, int col,
                char transa='N', char transb='N', int num=1,
                double beta=0, double alpha=1);
  
  /* SVD between first num and the rest indexes
         A = U * S * V
     returns the number of singular values (length of S) */
  int svd(tensor& U, tensor<double>& S, tensor& V, int num=1);

  int eig(double* val, T* vec);

};

#endif

#ifndef _MY_TENSOR_
#define _MY_TENSOR_

#include "global.h"

// Attention: the constructor function and functions:
// 'exchange', 'transpose', 'cut'
// can only accept at most 8 indexes

// T = 'double', 'complex<double>'

template <typename T> 
class tensor
{
  template <typename> friend class tensor;
  template <typename> friend class qtensor;

private:
  vector<int> index_;
  vector< T > val_;
  
  inline T* begin() { return &*val_.begin(); }
  
  int dim_() {
    if(index_.size()==0) return 0;
    int ret = 1;
    for(int i=0;i<index_.size();i++) ret *= index_[i];
    return ret;
  }
	
public:
  // claim space with all elements as zero
  tensor(int i0=0, int i1=0, int i2=0, int i3=0, 
         int i4=0, int i5=0, int i6=0, int i7=0);
		 
  // creat tensor with input values
  tensor(T* val,
         int i0, int i1=0, int i2=0, int i3=0, 
         int i4=0, int i5=0, int i6=0, int i7=0);

  void update(T val,
              int i0, int i1=-1, int i2=-1, int i3=-1, 
              int i4=-1, int i5=-1, int i6=-1, int i7=-1);

  // assign random values to the tensor
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

  tensor conjugate();

  // exchange the two indexes at positions a, b (start from 0)
  // for a two-tensor and a=0, b=1, this is just matrix transpose.
  tensor exchange(int a=0, int b=1) const;

  // transpose as a matrix with the first num indexes as left index
  // and the rest indexes as right index
  tensor shift(int num=1) const;
  
  // create a num-tensor with only diagonal elements non-zero
  tensor& diag(int num=1);
  
  // combine the index from min to max (including max)
  tensor& combine(int min=0, int max=1);
  
  // retain 'cut' number of index 'index'
  tensor resize(int index, int cutoff) const;
  
  // result = alpha * A + beta * B
  tensor& plus(const tensor& A, const tensor& B, T alpha=1, T beta=1);

  tensor times(T alpha) const {
    tensor<T> ret = *this;
    for(int i=0;i<val_.size();i++) ret.val_[i] *= alpha;
    return ret;
  }
  
  tensor operator+(const tensor& A) const {
    tensor<T> ret;
    ret.plus(*this, A);
    return ret;
  }
  
  // contract tensor A and B with the num left/right most index
  // the result is added to the original tensor (if not empty);
  tensor& contract(tensor& A, tensor& B, char transa='N', char transb='N',
                   int num=1, bool cover=false);
  
  tensor operator*(tensor& A) {
    tensor<T> ret;
    ret.contract(*this, A);
    return ret;
  }
  
  // contract index a of tensor A with index b of tensor B
  // the final index has the form
  // A0 A1 A2 ... Aa-1 B0 B1 ... Bm Aa+1 ... An
  // if cover = false, the result tensor will be added to the original one
  // if cover = true, the original tensor will be destroyed.
  tensor& contract(const tensor& A, int a, const tensor& B, int b,
                   bool cover=false, double sign = 1);
  
  // multiply the tensor by an array (as a matrix)
  // only multiply the 'num' left/right most index
  // for the array, use the same notation as matrix multiplication
  // beta=0 will erase the original data in out
  // bata=1 will add the result to the original data in out
  void contract(T* out, T* in, int row, int col,
                char transa='N', char transb='N', int num=1,
                double beta=0, double alpha=1);
  
  // SVD between first num and the rest indexes
  // A = U * S * V
  // returns the number of singular values (length of S)
  int svd(tensor& U, tensor<double>& S, tensor& V, int num=1);

  int eig(double* val, T* vec);

};

#endif

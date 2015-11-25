#ifndef _TENSOR_QUANTUM_
#define _TENSOR_QUANTUM_

#include "global.h"
#include "tensor_dense.h"
#include "tensor_sparse.h"

// Attention: the constructor function and functions:
// 'exchange', 'transpose', 'cutoff'
// can only accept at most 8 indexes


// ***********************************************************************
// the tensor with qunatum numbers, T = 'dtensor', 'stensor'
// can only accept at most 8 indexes

template <typename T> 
class qtensor
{
private:
  int sign_;
  vector<int> index_;  // the symmetry sector numbers for all indexes
  vector< vector<int> > sym_;  // we assume the symmetries are represented by 0,1,2...
  vector< T > val_;
	
public:
  // claim space with all elements as zero
  qtensor(int i0=1, int i1=0, int i2=0, int i3=0, 
         int i4=0, int i5=0, int i6=0, int i7=0);

  void update(T val,
              int i0, int i1=-1, int i2=-1, int i3=-1, 
              int i4=-1, int i5=-1, int i6=-1, int i7=-1);

  void print() const;

  // exchange the two indexes at positions a, b (start from 0)
  // for a two-tensor and a=0, b=1, this is just matrix transpose.
  qtensor exchange(int a=0, int b=1) const;

  // conjugate
  qtensor conjugate() const;

  // transpose as a matrix with the first num indexes as left index
  // and the rest indexes as right index
  qtensor transpose(int num=1) const;
  
  // result = alpha * A + beta * B
//  stensor& plus(const stensor& A, const stensor& B, T alpha=1, T beta=1);

  // create a num-tensor with only diagonal elements non-zero
  qtensor& diag(int num=1);
};
#endif

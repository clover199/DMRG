
#include "tensor_sparse.h"

template <typename T>
stensor<T>::stensor(int i0, int i1, int i2, int i3,
                  int i4, int i5, int i6, int i7)
{
  index_.clear();
  if(i0>0){ index_.push_back(i0); 
  if(i1>0){ index_.push_back(i1);
  if(i2>0){ index_.push_back(i2);
  if(i3>0){ index_.push_back(i3);
  if(i4>0){ index_.push_back(i4);
  if(i5>0){ index_.push_back(i5);
  if(i6>0){ index_.push_back(i6);
  if(i7>0){ index_.push_back(i7); }}}}}}}}
  else cout << "Error in tensor_sparse constructor: "
               "The input for the first index is less than 1.\n";
}


template <typename T>
void stensor<T>::update(T val,
                  int i0, int i1, int i2, int i3,
                  int i4, int i5, int i6, int i7)
{
  vector<int> loc;
  if(i0>=0){ loc.push_back(i0);
  if(i1>=0){ loc.push_back(i1);
  if(i2>=0){ loc.push_back(i2);
  if(i3>=0){ loc.push_back(i3);
  if(i4>=0){ loc.push_back(i4);
  if(i5>=0){ loc.push_back(i5);
  if(i6>=0){ loc.push_back(i6);
  if(i7>=0){ loc.push_back(i7);}}}}}}}}
  else cout << "Error in tensor_sparse update: "
               "The input for the first index is less than 0.\n";

  if( loc.size() == index_.size() )
  {
    cor_.push_back(loc);
    val_.push_back(val);
  }
  else
    cout << "Error in tensor_sparse update: "
            "The number of index doesn't match \n";
}


template <typename T>
int stensor<T>::a2n(const vector<int> in) const
{
  int result = 0;
  for(int i=0;i<in.size();i++) result = result*index_[i]+in[i];
  return result;
}


template <typename T>
vector<int> stensor<T>::n2a(int num) const
{
  vector<int> result(index_.size(),0);
  int i = result.size()-1;
  while(num && i>=0)
  {
    result[i] = num%index_[i];
    num = num/index_[i];
    i--;
  }
  return result;
}


template <typename T>
stensor<T>& stensor<T>::convertfrom(const dtensor<T>& in)
{
  index_ = in.index_;
  cor_.clear();
  val_.clear();
  for(int i=0;i<in.val_.size();i++) if(abs(in.val_[i])>ZERO)
  {
    cor_.push_back( n2a(i) );
    val_.push_back(in.val_[i]);
  }
  return *this;
}


template <typename T>
dtensor<T>& stensor<T>::convertto(dtensor<T>& in) const
{
  in.index_ = index_;
  in.val_.resize( dimension() );
  for(int i=0;i<in.val_.size();i++) in.val_[i] = 0;
  for(int i=0;i<val_.size();i++) in.val_[ a2n( cor_[i] ) ] += val_[i];
  return in;
}


template <typename T>
void stensor<T>::sort()
{
  dtensor<T> store;
  convertto(store);
  convertfrom(store);
}


template <typename T>
int stensor<T>::dimension() const 
{
  int dim = 1;
  for(int i=0;i<index_.size();i++) dim *= index_[i];
  return dim;
}


template <typename T>
int stensor<T>::dimension(int i) const { return index_[i]; }


template<typename T>
void stensor<T>::print() const
{
  cout << "Print all non-zero values of the tensor_sparse:\n"
       << "index=(" << index_[0];
  for(int i=1;i<index_.size();i++) cout << ", " << index_[i];
  cout << ") " << dimension() << "\n";
  for(int i=0;i<val_.size();i++) if(abs(val_[i])>ZERO)
  {
    cout << "(" << cor_[i][0];
    for(int j=1;j<index_.size();j++) cout << ", " << cor_[i][j];
    cout << ") " << val_[i] << "\n";
  }
}


template <typename T>
stensor<T> stensor<T>::exchange(int a, int b) const
{
  if( a>=index_.size() or b>=index_.size() or a<0 or b<0)
  {
    cout << "Error in tensor_sparse exchange: The input index number "
         << a << " and " << b << " is out of range 0~" << index_.size()-1 << endl;
    return *this;
  }
  
  stensor<T> ret(*this);
  ret.index_[a] = index_[b];
  ret.index_[b] = index_[a];
  
  for(int i=0;i<val_.size();i++)
  {
    int store = cor_[i][a];
    ret.cor_[i][a] = cor_[i][b];
    ret.cor_[i][b] = store;
  }
  ret.sort();
  return ret;
}


template<>
stensor<double> stensor<double>::conjugate() const
{
  return *this;
}


template<>
stensor< complex<double> > stensor< complex<double> >::conjugate() const
{
  stensor< complex<double> > ret;
  ret.index_ = index_;
  ret.cor_ = cor_;
  ret.val_.resize(val_.size());
  for(int i=0;i<val_.size();i++) ret.val_[i] = std::conj(val_[i]);
  return ret;
}


template <typename T>
stensor<T> stensor<T>::transpose(int num) const
{
  if( num>=index_.size() or num<=0)
  {
    cout << "Error in tensor_sparse transpose: The input index number "
         << num << " is out of range 1~" << index_.size()-1 << endl;
    return *this;
  }
  
  stensor<T> ret(*this);
  for(int i=0;i<index_.size();i++) ret.index_[i] = index_[(i+num)%index_.size()];
  for(int i=0;i<val_.size();i++) for(int j=0;j<index_.size();j++)
    ret.cor_[i][j] = cor_[i][(j+num)%index_.size()];
  ret.sort();
  return ret;
}


template <typename T>
stensor<T> stensor<T>::cut(int index, int cutoff) const
{
  if( index>=index_.size() or index<0)
  {
    cout << "Error in tensor_sparse cutoff: The input index number "
         << index << " is out of range 0~" << index_.size()-1 << endl;
    return *this;
  }
  
  if(cutoff>index_[index])
  {
    cout << "Warning in tensor_dense cutoff: the cut number "
         << cutoff << "is larger than the original number " << index_[index] << endl;
    return *this;
  }
  
  stensor<T> ret;
  ret.index_ = index_;
  ret.cor_.clear();
  ret.val_.clear();
  ret.index_[index] = cutoff;
  for(int i=0;i<val_.size();i++) if(cor_[i][index]<cutoff)
  {
    ret.cor_.push_back( cor_[i] );
    ret.val_.push_back( val_[i] );
  }
  return ret;
}


template <typename T>
stensor<T>& stensor<T>::plus(const stensor& A, const stensor& B, T alpha, T beta)
{
  if( A.index_.size() - B.index_.size() )
  {
    cout << "Error in tensor_sparse plus: "
            "Number of indexes doesn't match.\n" ;
    return *this;
  }

  int diff = 0;
  for(int i=0;i<A.index_.size();i++) diff += (A.index_[i]-B.index_[i]);
  if(diff)
  {
    cout << "Error in tensor_sparse plus: "
            "The indexes are not the same.\n";
    return *this;
  }
  
  index_ = A.index_;
  cor_ = A.cor_;
  val_ = A.val_;
  cor_.resize(A.cor_.size()+B.cor_.size());
  val_.resize(A.val_.size()+B.val_.size());
  for(int i=0;i<A.val_.size();i++) val_[i] *= alpha;
  for(int i=A.val_.size();i<val_.size();i++)
  {
    cor_[i] = B.cor_[i];
    val_[i] = beta * B.val_[i];
  }
  sort();
  return *this;
}


template <typename T>
stensor<T>& stensor<T>::diag(int num)
{
  int diff = 0;
  for(int i=1;i<index_.size();i++) diff += abs(index_[i]-index_[0]);
  if(diff)
  {
    cout << "Error in tensor_sparse diag: "
            "The indexes are not the same.\n";
    return *this;
  }
  
  vector<T> vals;  // used to store the diagonal elements
  vector<int> cor;  // used to store the corresponding coordinate
  vals.clear();
  cor.clear();
  for(int i=0;i<val_.size();i++)
  {
    diff = 0;
    for(int j=1;j<index_.size();j++) diff = abs(cor_[i][j]-cor_[i][0]);
    if(diff==0)
    {
      vals.push_back(val_[i]);
      cor.push_back(cor_[i][0]);
    }
  }
  index_.resize(num);
  for(int i=0;i<num;i++) index_[i] = index_[0];
  val_ = vals;
  cor_.clear();
  for(int i=0;i<vals.size();i++)
  {
    vector<int> temp; temp.clear();
    for(int j=0;j<num;j++) temp.push_back(cor[i]);
    cor_.push_back(temp);
  }
  return *this;
}


template class stensor<double>;
template class stensor< complex<double> >;


#include "tensor_dense.h"

template <typename T>
dtensor<T>::dtensor(int i0, int i1, int i2, int i3,
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
  else cerr << "Error in tensor_dense constructor: "
               "The input for the first index is less than 1.\n";
					 
  int dim = i0;
  for(int i=1;i<index_.size();i++) dim *= index_[i];
  val_.resize(dim);
  for(int i=0;i<dim;i++) val_[i] = 0;
}


template <typename T>
void dtensor<T>::update(T val,
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
  else cerr << "Error in tensor_dense update: "
               "The input for the first index is less than 0.\n";

  if( loc.size() == index_.size() )
  {
    int j = loc[0];
    for(int i=1;i<loc.size();i++) j = j*index_[i]+loc[i];
	val_[j] += val;
  }
  else
    cerr << "Error in tensor_dense update: "
            "The number of index doesn't match \n";
}


template<typename T>
int dtensor<T>::dimension() const { return val_.size(); }


template<typename T>
int dtensor<T>::dimension(int i) const { return index_[i]; }


template<typename T>
void dtensor<T>::print() const
{ 
  cout << "Print all non-zero values of the tensor_dense:\n"
       << "index=(" << index_[0];
  for(int i=1;i<index_.size();i++) cout << ", " << index_[i];
  cout << ") size=" << val_.size() << "\n";
  for(int i=0;i<val_.size();i++) if(abs(val_[i])>ZERO)
  {
    int store = i;
    int k = index_.size()-1;
    vector<int> temp(index_.size(),0);  // used to store the index
    while(store && k>=0)
    {
      temp[k] = store%index_[k];
      store = store/index_[k];
      k--;
    }
    cout << "(" << temp[0];
    for(int j=1;j<temp.size();j++) cout << ", " << temp[j];
    cout << ") " << val_[i] << "\n";
  }
}


template<typename T>
dtensor<T> dtensor<T>::exchange(int a, int b) const
{
  if( a>=index_.size() or b>=index_.size() or a<0 or b<0)
  {
    cerr << "Error in tensor_dense exchange: The input index number "
         << a << " and " << b << " is out of range 0~" << index_.size()-1 << endl;
    return *this;
  }
  
  dtensor<T> ret;
  ret.index_ = index_;
  ret.index_[a] = index_[b];
  ret.index_[b] = index_[a];

  vector<int> TOP(8,1);
  for(int i=0;i<index_.size();i++) TOP[i] = index_[i];
  ret.val_.resize(val_.size());
  vector<int> IND(8,0);
  for(IND[0]=0;IND[0]<TOP[0];IND[0]++)
  for(IND[1]=0;IND[1]<TOP[1];IND[1]++)
  for(IND[2]=0;IND[2]<TOP[2];IND[2]++)
  for(IND[3]=0;IND[3]<TOP[3];IND[3]++)
  for(IND[4]=0;IND[4]<TOP[4];IND[4]++)
  for(IND[5]=0;IND[5]<TOP[5];IND[5]++)
  for(IND[6]=0;IND[6]<TOP[6];IND[6]++)
  for(IND[7]=0;IND[7]<TOP[7];IND[7]++)
  {
    int loc1 = IND[0];
    for(int i=1;i<index_.size();i++) loc1 = loc1*index_[i]+IND[i];
    vector<int> index(IND);
    index[a] = IND[b];
    index[b] = IND[a];
    int loc2 = index[0];
    for(int i=1;i<index_.size();i++) loc2 = loc2*ret.index_[i]+index[i];
    ret.val_[loc2] = val_[loc1];
  }
  return ret;
}


template<>
dtensor<double> dtensor<double>::conjugate() const
{
  return *this;
}


template<>
dtensor< complex<double> > dtensor< complex<double> >::conjugate() const
{
  dtensor< complex<double> > ret;
  ret.index_ = index_;
  ret.val_.resize(val_.size());
  for(int i=0;i<val_.size();i++) ret.val_[i] = std::conj(val_[i]);
  return ret;
}


template<typename T>
dtensor<T> dtensor<T>::transpose(int num) const
{
  if( num>=index_.size() or num<=0)
  {
    cerr << "Error in tensor_dense transpose: The input index number "
         << num << " is out of range 1~" << index_.size()-1 << endl;
    return *this;
  }
  
  dtensor<T> ret;
  ret.index_.resize(index_.size());
  for(int i=0;i<index_.size();i++) ret.index_[i] = index_[(i+num)%index_.size()];

  vector<int> TOP(8,1);
  for(int i=0;i<index_.size();i++) TOP[i] = index_[i];
  ret.val_.resize(val_.size());
  vector<int> IND(8,0);
  for(IND[0]=0;IND[0]<TOP[0];IND[0]++)
  for(IND[1]=0;IND[1]<TOP[1];IND[1]++)
  for(IND[2]=0;IND[2]<TOP[2];IND[2]++)
  for(IND[3]=0;IND[3]<TOP[3];IND[3]++)
  for(IND[4]=0;IND[4]<TOP[4];IND[4]++)
  for(IND[5]=0;IND[5]<TOP[5];IND[5]++)
  for(IND[6]=0;IND[6]<TOP[6];IND[6]++)
  for(IND[7]=0;IND[7]<TOP[7];IND[7]++)
  {
    int loc1 = IND[0];
    for(int i=1;i<index_.size();i++) loc1 = loc1*index_[i]+IND[i];
    vector<int> index(IND);
    for(int i=0;i<index_.size();i++) index[i] = IND[(i+num)%index_.size()];
    int loc2 = index[0];
    for(int i=1;i<index_.size();i++) loc2 = loc2*ret.index_[i]+index[i];
    ret.val_[loc2] = val_[loc1];
  }
  return ret;
}


template <typename T>
dtensor<T> dtensor<T>::cut(int index, int cutoff) const
{
  if( index>=index_.size() or index<0)
  {
    cerr << "Error in tensor_dense cutoff: The input index number "
         << index << " is out of range 0~" << index_.size()-1 << endl;
    return *this;
  }
  
  if(cutoff>index_[index])
  {
    cerr << "Warning in tensor_dense cutoff: the cut number "
         << cutoff << "is larger than the original number " << index_[index] << endl;
    return *this;
  }
  
  dtensor<T> ret;
  ret.index_ = index_;
  ret.index_[index] = cutoff;

  vector<int> TOP(8,1);
  for(int i=0;i<index_.size();i++) TOP[i] = ret.index_[i];
  ret.val_.resize( val_.size() * ret.index_[index] / index_[index] );
  vector<int> IND(8,0);
  for(IND[0]=0;IND[0]<TOP[0];IND[0]++)
  for(IND[1]=0;IND[1]<TOP[1];IND[1]++)
  for(IND[2]=0;IND[2]<TOP[2];IND[2]++)
  for(IND[3]=0;IND[3]<TOP[3];IND[3]++)
  for(IND[4]=0;IND[4]<TOP[4];IND[4]++)
  for(IND[5]=0;IND[5]<TOP[5];IND[5]++)
  for(IND[6]=0;IND[6]<TOP[6];IND[6]++)
  for(IND[7]=0;IND[7]<TOP[7];IND[7]++)
  {
    int loc1 = IND[0];
    for(int i=1;i<index_.size();i++) loc1 = loc1*index_[i]+IND[i];
    int loc2 = IND[0];
    for(int i=1;i<index_.size();i++) loc2 = loc2*ret.index_[i]+IND[i];
    ret.val_[loc2] = val_[loc1];
  }
  return ret;
}


template<typename T>
dtensor<T>& dtensor<T>::plus(const dtensor& A, const dtensor& B, T alpha, T beta)
{
  if( A.index_.size() - B.index_.size() )
  {
    cerr << "Error in tensor plus: "
            "Number of indexes doesn't match.\n" ;
    return *this;
  }

  int diff = 0;
  for(int i=0;i<A.index_.size();i++) diff += (A.index_[i]-B.index_[i]);
  if(diff)
  {
    cerr << "Error in tensor plus: "
            "The indexes are not the same.\n";
    return *this;
  }
  
  index_ = A.index_;
  val_.resize(A.val_.size());
  for(int i=0;i<A.val_.size();i++) val_[i] = alpha * A.val_[i] + beta * B.val_[i];
  return *this;
}


#include "basic.h"

template<typename T>
dtensor<T>& dtensor<T>::contract(dtensor& A, dtensor& B, int num)
{
  if( num>A.index_.size() or num>B.index_.size() )
  {
    cerr << "Error in tensor_dense contract: "
            "Too many indexes to be contracted.\n" ;
    return *this;
  }

  int diff = 0;
  for(int i=0;i<num;i++) diff += (A.index_[A.index_.size()-num+i]-B.index_[i]);
  if(diff)
  {
    cerr << "Error in tensor_dense contract: "
            "The indexes are not the same.\n";
    return *this;
  }
  
  index_.clear();
  for(int i=0;i<A.index_.size()-num;i++) index_.push_back(A.index_[i]);
  for(int i=num;i<B.index_.size();i++) index_.push_back(B.index_[i]);
  int m = 1;
  for(int i=0;i<A.index_.size()-num;i++) m *= A.index_[i];
  int n = 1;
  for(int i=num;i<B.index_.size();i++) n *= B.index_[i];
  int k = 1;
  for(int i=0;i<num;i++) k *= B.index_[i];
  val_.resize(m*n);
  tgemm('N', 'N', m, n, k, A.begin(), B.begin(), begin(), 1, 0);
  return *this;
}


template<typename T>
dtensor<T>& dtensor<T>::combine(dtensor& A, dtensor& B)
{
  index_.resize( A.index_.size() + B.index_.size() );
  for(int i=0;i<A.index_.size();i++) index_[i] = A.index_[i];
  for(int i=0;i<B.index_.size();i++) index_[A.index_.size()+i] = B.index_[i];
  
  val_.resize( A.val_.size() * B.val_.size() );
  for(int i=0;i<A.val_.size();i++) for(int j=0;j<B.val_.size();j++)
    val_[i*B.val_.size()+j] = A.val_[i] * B.val_[j];
  return *this;
}


template<typename T>
dtensor<T>& dtensor<T>::diag(int num)
{
  int diff = 0;
  for(int i=1;i<index_.size();i++) diff += abs(index_[i]-index_[0]);
  if(diff)
  {
    cerr << "Error in tensor_dense diag: "
            "The indexes are not the same.\n";
    return *this;
  }
  
  vector<T> vals;  // used to store the diagonal elements
  vals.resize(index_[0]);
  for(int i=0;i<index_[0];i++)
  {
    int loc = i;
    for(int j=1;j<index_.size();j++) loc = loc*index_[j] + i;
    vals[i] = val_[loc];
  }
  index_.resize(num);
  for(int i=0;i<num;i++) index_[i] = vals.size();
  val_.resize(pow(vals.size(),num));
  for(int i=0;i<val_.size();i++) val_[i] = 0;
  for(int i=0;i<vals.size();i++)
  {
    int loc = i;
    for(int j=1;j<num;j++) loc = loc*index_[j] + i;
    val_[loc] = vals[i];
  }
  return *this;
}


template <typename T>
int dtensor<T>::svd(dtensor& U, dtensor<double>& S, dtensor& V, int num)
{
  if( num>=index_.size() or num<=0 )
  {
    cerr << "Error in tensor_dense svd: "
            "The number of left indexes should be 1~"<< index_.size()-1 << endl;
    return 0;
  }
  
  int left = 1;
  for(int i=0;i<num;i++) left *= index_[i];
  int right = 1;
  for(int i=num;i<index_.size();i++) right *= index_[i];
  int ret = left < right ? left : right;  // min{left,right}
  
  U.index_.resize(num+1);
  for(int i=0;i<num;i++) U.index_[i] = index_[i];
  U.index_[num] = left;
  U.val_.resize(left*left);
  V.index_.resize(index_.size()-num+1);
  for(int i=0;i<index_.size()-num;i++) V.index_[i+1] = index_[i+num];
  V.index_[0] = right;
  V.val_.resize(right*right);
  S.index_.resize(1);
  S.index_[0] = ret;
  S.val_.resize(ret);
  
  tgesvd(begin(), left, right, S.begin(), U.begin(), V.begin());
  return ret;
}
  
/* Explicit instantiation */
template class dtensor<double>;
template class dtensor< complex<double> >;

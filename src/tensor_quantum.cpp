
#include "tensor_quantum.h"


template <typename T>
qtensor<T>::qtensor(int i0, int i1, int i2, int i3,
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
  else cerr << "Error in tensor_quantum constructor: "
               "The input for the first index is less than 1.\n";
}


template <typename T>
void qtensor<T>::update(T val,
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
  else cerr << "Error in tensor_quantum update: "
               "The input for the first index is less than 0.\n";

  if( loc.size() == index_.size() )
  {
    sym_.push_back(loc);
    val_.push_back(val);
  }
  else
    cerr << "Error in tensor_sparse update: "
            "The number of index doesn't match \n";
}


template<typename T>
void qtensor<T>::print() const
{
  cout << "Print all non-zero symmetry sectors of the tensor_quantum:\n"
       << "index=(" << index_[0];
  for(int i=1;i<index_.size();i++) cout << ", " << index_[i];
  cout << ")\n";
  for(int i=0;i<val_.size();i++)
  {
    cout << "i (" << sym_[i][0];
    for(int j=1;j<index_.size();j++) cout << ", " << sym_[i][j];
    cout << ") dimension:(" << val_[i].dimension(0);
    for(int j=1;j<val_[i].rank();j++) cout << ". " << val_[i].dimension(j);
    cout << ")\n";
  }
}


template <typename T>
qtensor<T> qtensor<T>::exchange(int a, int b) const
{
  if( a>=index_.size() or b>=index_.size() or a<0 or b<0)
  {
    cerr << "Error in tensor_quantum exchange: The input index number "
         << a << " and " << b << " is out of range 0~" << index_.size()-1 << endl;
    return *this;
  }
  
  qtensor<T> ret;  
  ret.sign_ = sign_;
  ret.index_ = index_;
  ret.sym_ = sym_;
  ret.val_.resize(val_.size());
  for(int i=0;i<val_.size();i++) ret.val_[i] = val_[i].exchange(a,b);
  return ret;
}


template<typename T>
qtensor<T> qtensor<T>::conjugate() const
{
  qtensor<T> ret;
  ret.sign_ = sign_;
  ret.index_ = index_;
  ret.sym_ = sym_;
  ret.val_.resize(val_.size());
  for(int i=0;i<val_.size();i++) ret.val_[i] = val_[i].conjugate();
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::transpose(int num) const
{
  if( num>=index_.size() or num<0)
  {
    cerr << "Error in tensor_quantum transpose: The input index number "
         << num << " is out of range 1~" << index_.size()-1 << endl;
    return *this;
  }
  
  qtensor<T> ret;
  ret.sign_ = sign_;
  ret.index_ = index_;
  ret.sym_ = sym_;
  ret.val_.resize(val_.size());
  for(int i=0;i<val_.size();i++) ret.val_[i] = val_[i].transpose(num);
  return ret;
}


template <typename T>
qtensor<T>& qtensor<T>::diag(int num)
{
  
}

template class qtensor< dtensor<double> >;
template class qtensor< dtensor< complex<double> > >;
template class qtensor< stensor<double> >;
template class qtensor< stensor< complex<double> > >;

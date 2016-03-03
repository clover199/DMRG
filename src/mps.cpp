#include "mps.h"


template <typename T> 
void mps<T>::print() const 
{
  cout << "Print all MPS state: (cut=" << cut_ << ")\n";
  for(int i=0;i<state_.size();i++)
  {
    cout << i << ": ";cout<<endl;
    state_[i].print();
    cout << endl;
  }
  cout << "Print all stored tensor" << endl;
  for(int i=0;i<store_.size();i++)
  {
    cout << i << ": ";
    store_[i].print();
    cout << endl;
  }
#ifdef PBC
  cout << "Print all edge tensor" << endl;
  for(int i=0;i<edge_.size();i++)
  {
    cout << i << ": ";
    edge_[i].print();
    cout << endl;
  }
#endif
}

template <>
void mps<double>::prep_calc()
{
  store_.clear();
#ifdef PBC
  edge_.clear();
#endif
  state_[cut_] = singular_ * state_[cut_];
}


template <>
void mps<complex<double> >::prep_calc()
{
  if(store_.size()) state_[cut_] = singular_.comp() * state_[cut_];
  else cout << "Warning in mps prep_calc: already done.\n";
  store_.clear();
#ifdef PBC
  edge_.clear();
#endif
}


template <typename T>
mps<T> mps<T>::inverse() const 
{
  mps<T> ret;
  int size = state_.size();
  ret.cut_ = size-cut_;
  ret.state_.resize(size);
  for(int i=1;i<size-1;i++) ret.state_[size-1-i] = state_[i].exchange(0,2);
  ret.state_[0] = state_[size-1].exchange(0,1);
  ret.state_[size-1] = state_[0].exchange(0,1);
  ret.singular_ = singular_;
  return ret;
}


template <>
void mps<double>::move_cut_left()
{
  if(cut_<2)
  {
    cerr << "Error in mps move_cut_left: cut=" << cut_ << endl;
    return;
  }
  qtensor<double> temp;
  temp = state_[cut_-2] * state_[cut_-1] * singular_;
  if(cut_>2) temp.svd( state_[cut_-2], singular_, state_[cut_-1], 2);
  else temp.svd( state_[cut_-2], singular_, state_[cut_-1], 1);
  cut_ --;
}


template <>
void mps<complex<double> >::move_cut_left()
{
  if(cut_<2)
  {
    cerr << "Error in mps move_cut_left: cut=" << cut_ << endl;
    return;
  }
  qtensor<complex<double> > temp, store;
  store = singular_.comp();
  temp = state_[cut_-2] * state_[cut_-1] * store;
  if(cut_>2) temp.svd( state_[cut_-2], singular_, state_[cut_-1], 2);
  else temp.svd( state_[cut_-2], singular_, state_[cut_-1], 1);
  cut_ --;
}

template <>
void mps<double>::move_cut_right()
{
  if(cut_>state_.size()-2)
  {
    cerr << "Error in mps move_cut_right: cut=" << cut_ << endl;
    return;
  }
  qtensor<double> temp;
  temp = singular_ * state_[cut_] * state_[cut_+1];
  temp.svd( state_[cut_], singular_, state_[cut_+1], 2);
  cut_ ++;
}

template <>
void mps<complex<double> >::move_cut_right()
{
  if(cut_>state_.size()-2)
  {
    cerr << "Error in mps move_cut_right: cut=" << cut_ << endl;
    return;
  }
  qtensor<complex<double> > temp;
  temp = singular_.comp() * state_[cut_] * state_[cut_+1];
  temp.svd( state_[cut_], singular_, state_[cut_+1], 2);
  cut_ ++;
}

template class mps<double>;
template class mps<complex<double> >;

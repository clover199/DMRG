
#include "tensor.h"
#include "basic.h"
#include "useful.h"

template <typename T>
tensor<T>::tensor(int i0, int i1, int i2, int i3,
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

  int dim = dim_();
  val_.resize(dim);
  for(int i=0;i<dim;i++) val_[i] = 0;
}


template <typename T>
tensor<T>::tensor(T* val,
                  int i0, int i1, int i2, int i3,
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
  else cerr << "Error in tensor constructor: "
               "The input for the first index is < 1.\n";

  int dim = dim_();
  val_.resize(dim);
  for(int i=0;i<dim;i++) val_[i] = val[i];
}


template <typename T>
void tensor<T>::update(T val,
                       int i0, int i1, int i2, int i3,
                       int i4, int i5, int i6, int i7)
{
  vector<int> loc; loc.clear();
  if(i0>=0){ loc.push_back(i0);
  if(i1>=0){ loc.push_back(i1);
  if(i2>=0){ loc.push_back(i2);
  if(i3>=0){ loc.push_back(i3);
  if(i4>=0){ loc.push_back(i4);
  if(i5>=0){ loc.push_back(i5);
  if(i6>=0){ loc.push_back(i6);
  if(i7>=0){ loc.push_back(i7);}}}}}}}}
  
  if( loc.size() == index_.size() )
  {
    bool check = true;
    for(int i=0;i<loc.size();i++) if(loc[i]>=index_[i])
    {
      cerr << "Error in tensor update: "
              "The input " << i+1 << "th index " << loc[i] 
           << " is out of range 0~" << index_[i]-1 << endl;
      check = false;
    }
    if(check)
    {
      int j = 0;
      for(int i=0;i<loc.size();i++) j = j*index_[i]+loc[i];
	val_[j] = val;
    }
  }
  else
    cerr << "Error in tensor update: The number of index doesn't match \n";
}


template <>
void tensor<double>::rand()
{
  for(int i=0;i<val_.size();i++) val_[i] = 10.0*std::rand()/RAND_MAX-5;
}


template <>
void tensor< complex<double> >::rand()
{
  for(int i=0;i<val_.size();i++)
    val_[i] = complex<double> (10.0*std::rand()/RAND_MAX-5, 
                               10.0*std::rand()/RAND_MAX-5);
}


template <typename T>
void tensor<T>::print(bool all) const
{
  cout.precision(10);
  vector< vector<int> > map;
  generate_map(map, index_);
  cout << "Print all non-zero values of the tensor:\n"
       << "dimension=(" << index_[0];
  for(int i=1;i<index_.size();i++) cout << ", " << index_[i];
  cout << ") size=" << val_.size() << "\n";

  for(int i=0;i<val_.size();i++) if(all)
  {
    cout << "(" << map[i][0];
    for(int j=1;j<index_.size();j++) cout << ", " << map[i][j];
    cout << ") " << val_[i] << "\n";
  }
  else if(abs(val_[i])>TOL)
  {
    cout << "(" << map[i][0];
    for(int j=1;j<index_.size();j++) cout << ", " << map[i][j];
    cout << ") " << val_[i] << "\n";
  }
}


template <typename T>
void tensor<T>::print_matrix() const
{
  cout.precision(3);
  cout << "Print the tensor as a matrix:\n";
  if(2==index_.size())
  {
    cout << "\t row=" << index_[0] << "  col=" << index_[1] << endl;
    for(int i=0;i<index_[0];i++)
    {
      for(int j=0;j<index_[1];j++)
      {
        T val = val_[i*index_[1]+j];
        if(abs(val)>ZERO) cout << val << "\t";
        else cout << (T)0 << "\t";
      }
      cout << endl;
    }
  }
  else
    cout << "This " << index_.size() << "-tensor cannot "
            "be printed as a matrix\n";
}


template <>
tensor<double> tensor<double>::conjugate()
{
  return *this;
}


template <>
tensor< complex<double> > tensor< complex<double> >::conjugate()
{
  for(int i=0;i<val_.size();i++) val_[i] = std::conj(val_[i]);
  return *this;
}


template <typename T>
tensor<T> tensor<T>::exchange(int a, int b) const
{
  if( a>=index_.size() or b>=index_.size() or a<0 or b<0)
  {
    cerr << "Error in tensor exchange: The input index number "<< a 
         << " and/or " << b << " is out of range 0~" << index_.size()-1
         << endl;
    return *this;
  }
  
  tensor<T> ret;
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


template <typename T>
tensor<T> tensor<T>::shift(int num) const
{
  if( num>=index_.size() or num<=0)
  {
    cerr << "Error in tensor transpose: The input index number "
         << num << " is out of range 1~" << index_.size()-1 << endl;
    return *this;
  }
  
  tensor<T> ret;
  ret.index_.resize(index_.size());
  for(int i=0;i<index_.size();i++)
    ret.index_[i] = index_[(i+num)%index_.size()];

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
tensor<T>& tensor<T>::diag(int num)
{
  int min = index_[0];
  for(int i=0;i<index_.size();i++) if(index_[i]<min) min = index_[i];
  
  vector<T> vals;  // used to store the diagonal elements
  vals.resize(min);
  for(int i=0;i<min;i++)
  {
    int loc = i;
    for(int j=1;j<index_.size();j++) loc = loc*index_[j] + i;
    vals[i] = val_[loc];
  }
  index_.resize(num);
  for(int i=0;i<num;i++) index_[i] = min;
  val_.resize(dim_());
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
tensor<T>& tensor<T>::combine(int min, int max)
{
  if(min==max) return *this;
  if( 0<=min and min<max and max<index_.size() )
  {
    vector<int> store(index_);
    index_.resize( store.size()-(max-min) );
    for(int i=min+1;i<=max;i++) index_[min] *= store[i];
    for(int i=1;i<store.size()-max;i++) index_[min+i] = store[max+i];
  }
  else
    cerr << "Error in tensor combine: "
            "The inputs do not satisfy 0 <= " << min << " < " << max
         << " <= " << index_.size()-1 << endl;
  return *this;
}


template <typename T>
tensor<T> tensor<T>::resize(int index, int cutoff) const
{
  if( index>=index_.size() or index<0)
  {
    cerr << "Error in tensor cutoff: The input index number "
         << index << " is out of range 0~" << index_.size()-1 << endl;
    return *this;
  }
  
  tensor<T> ret;
  ret.index_ = index_;
  ret.index_[index] = cutoff;

  vector<int> TOP(8,1);
  for(int i=0;i<index_.size();i++) TOP[i] = ret.index_[i];
  if(index_[index]<ret.index_[index]) TOP[index] = index_[index];
  ret.val_.resize( val_.size() * ret.index_[index] / index_[index] );
  for(int i=0;i<ret.val_.size();i++) ret.val_[i] = 0;
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


template <typename T>
tensor<T>& tensor<T>::plus(const tensor& A, const tensor& B,
                            T alpha, T beta)
{
  if( A.index_.size() - B.index_.size() )
  {
    cerr << "Error in tensor plus: Number of indexes doesn't match. "
         << A.index_.size() << "!=" << B.index_.size() << endl;
    return *this;
  }

  int diff = 0;
  for(int i=0;i<A.index_.size();i++) diff += (A.index_[i]-B.index_[i]);
  if(diff)
  {
    cerr << "Error in tensor plus: The indexes are not the same.\n";
    return *this;
  }
  
  index_ = A.index_;
  val_.resize(A.val_.size());
  for(int i=0;i<A.val_.size();i++)
    val_[i] = alpha * A.val_[i] + beta * B.val_[i];
  return *this;
}


template <typename T>
tensor<T>& tensor<T>::contract(tensor& A, tensor& B, char transa, char transb,
                               int num, bool cover)
{
  int indexa = 0;
  if( 'N'==transa or 'n'==transa ) indexa = A.index_.size()-num;
  int indexb = B.index_.size()-num;
  if( 'N'==transb or 'n'==transb ) indexb = 0;

  int diff = 0;
  for(int i=0;i<num;i++) diff += abs(A.index_[indexa+i]-B.index_[indexb+i]);
  if(diff)
  {
    cerr << "Error in tensor contract: The indexes are not the same.\n"
         << "A: ";
    for(int i=0;i<num;i++) cerr << A.index_[indexa+i] << " ";
    cerr << "\nB: ";
    for(int i=0;i<num;i++) cerr << B.index_[indexb+i] << " ";
    cerr << endl;
    return *this;
  }

  index_.resize(A.index_.size()+B.index_.size()-2*num);
  int indexa2 = num;
  if( 'N'==transa or 'n'==transa ) indexa2 = 0;
  int indexb2 = 0;
  if( 'N'==transb or 'n'==transb ) indexb2 = num;
  for(int i=0;i<A.index_.size()-num;i++)
    index_[i] = A.index_[indexa2+i];
  for(int i=0;i<B.index_.size()-num;i++)
    index_[A.index_.size()-num+i] = B.index_[indexb2+i];
  
  int k = 1;
  for(int i=0;i<num;i++) k *= A.index_[indexa+i];
  int m = A.val_.size()/k;
  int n = B.val_.size()/k;

  if(val_.size()==0) val_.resize(m*n);
  else if(cover) val_ = vector<T> (m*n, 0);
  else if(val_.size()!=m*n)
  {
    cerr << "Warning in tensor contract: val_ size doesn't match "
         << val_.size() << "!=" << m*n <<". Erase old data\n";
    val_ = vector<T> (m*n, 0);
  }
  zgemm(transa, transb, m, n, k, A.begin(), B.begin(), begin(), 1, 1);
  return *this;
}


template <typename T>
tensor<T>& tensor<T>::contract(const tensor& A, int a, const tensor& B, int b,
                               bool cover)
{
  if(A.index_[a]!=B.index_[b])
  {
    cerr << "Error in tensor contract: Dimensions do not match. "
         << A.index_[a] <<"!=" << B.index_[b] << endl;
    return *this;
  }
  int mid = A.index_[a];
  
  index_.resize(A.index_.size()+B.index_.size()-2);
  for(int i=0;i<a;i++) index_[i] = A.index_[i];
  for(int i=a;i<a+b;i++) index_[i] = B.index_[i-a];
  for(int i=a+b;i<a+B.index_.size()-1;i++)
    index_[i] = B.index_[i-a+1];
  for(int i=a+B.index_.size()-1;i<index_.size();i++)
    index_[i] = A.index_[i-B.index_.size()+2];

  int al = 1;
  for(int i=0;i<a;i++) al *= A.index_[i];
  int ar = 1;
  for(int i=a+1;i<A.index_.size();i++) ar *= A.index_[i];
  int bl = 1;
  for(int i=0;i<b;i++) bl *= B.index_[i];
  int br = 1;
  for(int i=b+1;i<B.index_.size();i++) br *= B.index_[i];

  int dim = dim_();
  if(val_.size()==0) val_.resize(dim);
  else if(cover) val_ = vector<T> (dim, 0);
  else if(val_.size()!=dim)
  {
    cerr << "Warning in tensor contract: val_ size doesn't match "
         << val_.size() << "!=" << dim <<". Erase old data\n";
    val_ = vector<T> (dim, 0);
  }
  for(int k=0;k<mid;k++)
  for(int i1=0;i1<al;i1++) for(int j1=0;j1<bl;j1++)
  for(int i2=0;i2<ar;i2++) for(int j2=0;j2<br;j2++)
    val_[i1*bl*br*ar+j1*br*ar+j2*ar+i2] += 
        A.val_[i1*mid*ar+k*ar+i2] * B.val_[j1*mid*br+k*br+j2];
  return *this;
}


template <typename T>
void tensor<T>::contract(T* out, T* in, int row, int col,
                         char transa, char transb, int num, double beta)
{
  int indexa = 0;
  if( 'N'==transa or 'n'==transa ) indexa = index_.size()-num;
 
  int m = val_.size();
  for(int i=0;i<num;i++) m /= index_[indexa+i];
  int n = row;
  int k = col;
  if( 'N'==transb or 'n'==transb ) n = col, k = row;
  
  if( (val_.size()/m) != k )
  {
    cerr << "Error in tensor contract: The indexes are not the same. "
         << val_.size()/m << "!=" << k << endl;
  }
  zgemm(transa, transb, m, n, k, begin(), in, out, 1, beta);
}


template <typename T>
int tensor<T>::svd(tensor& U, tensor<double>& S, tensor& V, int num)
{
  if( num>=index_.size() or num<=0 )
  {
    cerr << "Error in tensor svd: "
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

  bool check = true;  // if all elements are zero
  for(int i=0;i<val_.size();i++) if(abs(val_[i])>TOL) check = false;
  if(check)
  {
    for(int i=0;i<U.val_.size();i++) U.val_[i] = 0;
    for(int i=0;i<left;i++) U.val_[i*left+i] = 1;
    for(int i=0;i<V.val_.size();i++) V.val_[i] = 0;
    for(int i=0;i<right;i++) V.val_[i*right+i] = 1;
    for(int i=0;i<S.val_.size();i++) S.val_[i] = 0;
    return ret;
  }

  if(left*right==1)
  {
    U.val_[0] = 1;
    V.val_[0] = val_[0]/abs(val_[0]);
    S.val_[0] = abs(val_[0]);
    return 1;
  }
  zgesvd(begin(), left, right, S.begin(), U.begin(), V.begin());
  return ret;
}


template <typename T>
int tensor<T>::eig(double* val, T* vec)
{
  if(index_.size()!=2)
  {
    cerr << "Error in tensor eig: "
            "This " << index_.size() << "-tensor is not a matrix.\n";
    return 0;
  }
  if( index_[0]!=index_[1] )
  {
    cerr << "Error in tensor eig: This is not a square matrix. "
         << index_[0] << "!=" << index_[1] << endl;
    return 0;
  }
  int d = index_[0];
  zheev(begin(), d, val, vec);
  return d;
}


/* Explicit instantiation */
template class tensor<double>;
template class tensor< complex<double> >;

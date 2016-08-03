
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
tensor<T>::tensor(const vector<int>& index) 
{
  index_ = index;
  val_.resize(dim_());
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
               "The first index " << i0 << " should be non-zero.\n";

  int dim = dim_();
  val_.resize(dim);
  for(int i=0;i<dim;i++) val_[i] = val[i];
}
  
  
template <typename T>
tensor<T>::tensor(T* val, const vector<int>& index)
{
  index_ = index;
  int dim = dim_();
  val_.resize(dim);
  for(int i=0;i<dim;i++) val_[i] = val[i];
}


template <typename T>
void tensor<T>::update(T val,
                       int i0, int i1, int i2, int i3,
                       int i4, int i5, int i6, int i7)
{
  if(abs(val)<TOL) return;
  vector<int> loc;
  if(i0>=0){ loc.push_back(i0);
  if(i1>=0){ loc.push_back(i1);
  if(i2>=0){ loc.push_back(i2);
  if(i3>=0){ loc.push_back(i3);
  if(i4>=0){ loc.push_back(i4);
  if(i5>=0){ loc.push_back(i5);
  if(i6>=0){ loc.push_back(i6);
  if(i7>=0){ loc.push_back(i7);}}}}}}}}
  
  if( loc.size()==1)
  {
    val_.at(i0) += val;
    return ;
  }
  
  if( loc.size() != index_.size() )
  {
    cerr << "Error in tensor update: The number of input index " << loc.size() 
         << " doesn't match this " << index_.size() << "-tensor.\n";
    return ;
  }
  for(int i=0;i<loc.size();i++) if(loc[i]>=index_[i])
  {
    cerr << "Error in tensor update: The input " << i+1 << "th index " << loc[i] 
         << " is out of range 0~" << index_[i]-1 << endl;
    return;
  }
  
  val_[loc_(loc)] += val;
}


template <typename T>
void tensor<T>::zeros()
{
  for(int i=0;i<val_.size();i++) val_[i] = 0.0;
}


template <typename T>
void tensor<T>::ones()
{
  for(int i=0;i<val_.size();i++) val_[i] = 1.0;
}


template <>
void tensor<int>::rand()
{
  for(int i=0;i<val_.size();i++) val_[i] = std::rand()-RAND_MAX/2;
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
void tensor<T>::id(int a1, int b1, int a2, int b2,
                   int a3, int b3, int a4, int b4)
{
  vector<int> pair;
  if(a1>=0){ pair.push_back(a1);
  if(b1>=0){ pair.push_back(b1);
  if(a2>=0){ pair.push_back(a2);
  if(b2>=0){ pair.push_back(b2);
  if(a3>=0){ pair.push_back(a3);
  if(b3>=0){ pair.push_back(b3);
  if(a4>=0){ pair.push_back(a4);
  if(b4>=0){ pair.push_back(b4);}}}}}}}}

  for(int i=0;i<pair.size();i++) if(pair[i]>=index_.size())
  {
    cerr << "Error in tensor id: The input " << i+1 << "th index " << pair[i] 
         << " is out of range 0~" << index_.size()-1 << endl;
    return;
  }
  if(pair.size()%2)
  {
    cerr << "Error in tensor id: The input indexes should be couples "<< endl;
    return;
  }
  
  ones();
  vector< vector<int> > map;
  generate_map(map, index_);
  for(int i=0;i<val_.size();i++) for(int p=0;p<pair.size()/2;p++)
    if( map[i].at(pair[2*p])!=map[i].at(pair[2*p+1]) ) val_[i] = 0.0;
}


template <typename T>
bool tensor<T>::empty() const
{
  double sum = 0.0;
  for(int i=0;i<val_.size();i++) sum += abs(val_[i]);
  if(sum<TOL) return true;
  else return false;
}


template <typename T>
void tensor<T>::print() const
{
  cout.precision(10);
  cout << "Print all non-zero values of the tensor:\n"
       << "dimension=(";
  for(int i=0;i<index_.size();i++)
    if(i==0) cout << index_[i];
    else cout << ", " << index_[i];
  cout << ") size=" << val_.size() << "\n";

  vector<int> loc(index_.size(), 0);
  for(int i=0;i<val_.size();i++)
  {
    if(abs(val_[i])>TOL)
    {
      cout << "(";
      for(int j=0;j<index_.size();j++)
        if(j==0) cout << loc[j];
        else cout << ", " << loc[j];
      cout << ") " << val_[i] << "\n";
    }
    loc.at(loc.size()-1) += 1;
    refine_(loc);
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
        T val = val_.at(i*index_[1]+j);
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

template <typename T>
void tensor<T>::print_index() const
{
  cout << "(";
  for(int i=0;i<index_.size();i++)
    if(i==0) cout << index_[i];
    else cout << ", " << index_[i];
  cout << ")" << endl;
}



template <typename T>
void tensor<T>::write_file(fstream& name) const
{
  name.precision(32);
  if(!name.is_open()) cerr << "Error in write file: cannot open file\n";
  name << index_.size();
  for(int i=0;i<index_.size();i++) name << " " << index_[i];
  for(int i=0;i<dim_();i++) name << " " << val_[i];
}


template <typename T>
void tensor<T>::read_file(fstream& name)
{
  if(!name.is_open()) cerr << "Error in read file: cannot open file\n";
  int num;
  name >> num;
  index_.resize(num);
  for(int i=0;i<index_.size();i++) name >> index_[i];
  num = dim_();
  val_.resize(num);
  for(int i=0;i<val_.size();i++) name >> val_[i];
}


template <typename T>
void tensor<T>::write_file(ofstream& name) const
{
  name.precision(32);
  if(!name.is_open()) cerr << "Error in write file: cannot open file\n";
  name << index_.size();
  for(int i=0;i<index_.size();i++) name << " " << index_[i];
  for(int i=0;i<dim_();i++) name << " " << val_[i];
}


template <typename T>
void tensor<T>::read_file(ifstream& name)
{
  if(!name.is_open()) cerr << "Error in read file: cannot open file\n";
  int num;
  name >> num;
  index_.resize(num);
  for(int i=0;i<index_.size();i++) name >> index_[i];
  num = dim_();
  val_.resize(num);
  for(int i=0;i<val_.size();i++) name >> val_[i];
}


template <>
void tensor<double>::conjugate() { return ; }

template <>
void tensor< complex<double> >::conjugate()
{
  for(int i=0;i<val_.size();i++) val_[i] = std::conj(val_[i]);
}


template <>
void tensor<double>::to_double(tensor<double>& ret) const
{
  ret.index_ = index_;
  ret.val_ = val_;
}

template <>
void tensor<complex<double> >::to_double(tensor<double>& ret) const
{ 
  cerr << "Warning in tensor to_double: discard imaginary part. " << endl;
  ret.index_ = index_;
  ret.val_.resize(val_.size());
  for(int i=0;i<val_.size();i++) ret.val_[i] = val_[i].real();
}


template <typename T>
void tensor<T>::to_complex(tensor<complex<double> >& ret) const
{
  ret.index_ = index_;
  ret.val_.resize(val_.size());
  for(int i=0;i<val_.size();i++) ret.val_[i] = val_[i]*1.0;
}


template <typename T>
void tensor<T>::exchange(tensor& ret, int a, int b) const
{  
  if( a>=index_.size() or b>=index_.size() or a<0 or b<0)
  {
    cerr << "Error in tensor exchange: The input index number "<< a 
         << " and/or " << b << " is out of range 0~" << index_.size()-1 << endl;
    return ;
  }
  
  ret.index_ = index_;

  if(a==b)
  {
    ret.val_ = val_;
    return ;
  }
  
  ret.index_[a] = index_[b];
  ret.index_[b] = index_[a];
  
  ret.val_.resize(val_.size());
  vector< vector<int> > map;
  generate_map(map, ret.index_);
  for(int i=0;i<map.size();i++)
  {
    int temp = map[i][a];
    map[i][a] = map[i][b];
    map[i][b] = temp;
    ret.val_[i] = val_.at( loc_(map[i]) );
  }
}


template <typename T>
void tensor<T>::shift(tensor& ret, int num) const
{
  if( num>index_.size() or num<0)
  {
    cerr << "Error in tensor transpose: The input index number "
         << num << " is out of range 0~" << index_.size() << endl;
    return ;
  }
  
  if(num==0 or num==index_.size())
  {
    ret.index_ = index_;
    ret.val_ = val_;
    return ;
  }

  ret.index_.resize(index_.size());
  for(int i=0;i<index_.size();i++)
    ret.index_[i] = index_[(i+num)%index_.size()];

  ret.val_.resize(val_.size());
  vector< vector<int> > map;
  generate_map(map, ret.index_);
  vector<int> index(index_.size(),0);
  for(int i=0;i<map.size();i++)
  {
    for(int j=0;j<index.size();j++) index[j] = map[i][(j+num)%index_.size()];
    ret.val_[i] = val_.at( loc_(index) );
  }
}


template <typename T>
void tensor<T>::take(tensor& ret, int n, int k) const
{
  if(n<0 or n>=index_.size())
  {
    cout << "Error in tensor take: "
            "the index " << n << " is out of the range 0~"
         << index_.size()-1 << endl;
    return ;
  }
  if( k<0 or k>=index_[n])
  {
    cout << "Error in tensor take: "
            "the number " << k << " is out of the dimension range 0~" 
         << index_[n] << endl;
    return ;
  }
  
  ret.index_.resize(index_.size()-1);
  for(int i=0;i<n;i++) ret.index_[i] = index_[i];
  for(int i=n;i<ret.index_.size();i++) ret.index_[i] = index_[i+1];

  int r = 1;
  for(int i=n+1;i<index_.size();i++) r *= index_[i];
  ret.val_.resize(dimension()/dimension(n));
  for(int i=0;i<ret.val_.size();i++)
  {
    int j = ( (i/r) * index_[n] + k ) * r + (i%r) ;
    ret.val_.at(i) = val_.at(j);
  }
}


template <typename T>
void tensor<T>::diag(int num)
{
  int min = index_.at(0);
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
}


template <typename T>
void tensor<T>::combine(int min, int max)
{
  if(min==max) return ;
  
  if( min<0 or max<min or max>=index_.size() )
  {
    cerr << "Error in tensor combine: "
            "The inputs do not satisfy 0 <= " << min << " < " << max
         << " <= " << index_.size()-1 << endl;
    return ;
  }
  
  for(int i=min+1;i<=max;i++) index_[min] *= index_[i];
  for(int i=1;i<index_.size()-max;i++) index_[min+i] = index_[max+i];
  int new_size = index_.size() - (max-min);
  index_.resize(new_size);
}


template <typename T>
void tensor<T>::resize(int k, int n)
{
  if(index_.size()==0) return;
  
  if( k>=index_.size() or k<0)
  {
    cerr << "Error in tensor resize: The input index number "
         << k << " is out of range 0~" << index_.size()-1 << endl;
    return ;
  }

  if(n==index_[k]) return;

  int rk = 1;  // radix for indexes before k (before resize)
  for(int i=k+1;i<index_.size();i++) rk *= index_[i];
  int rk1 = rk * n; // radix for indexes before k (after resize)
  rk *= index_[k];

  int n_old = index_[k];
  index_[k] = n;
  
  vector<T> val(val_);
  val_.resize(dim_());
  for(int i=0;i<val_.size();i++)
  {
    int l = i/rk1;
    int r = i%rk1;
    if(r*n/rk1 < n_old) val_.at(i) = val.at(l*rk+r);
    else val_.at(i) = 0.0;
  }
}


template <typename T>
void tensor<T>::resize(const vector<int>& index)
{
  index_ = index;
  val_.resize(dim_());
  for(int i=0;i<val_.size();i++) val_[i] = 0.0;
}


template <typename T>
T tensor<T>::trace() const
{
  T val = 0.0;
  if(index_.size()==0) return val;
  
  int min = index_[0];
  for(int i=0;i<index_.size();i++) if(index_[i]<min) min = index_[i];
  
  for(int i=0;i<min;i++)
  {
    int loc = i;
    for(int j=1;j<index_.size();j++) loc = loc*index_[j] + i;
    val += val_.at(loc);
  }
  return val;
}


template <typename T>
T tensor<T>::trace(tensor& A) const
{
  if(index_!=A.index_)
  {
    cerr << "Error in tensor trace: The indexes are not the same.\n";
    return 0.0;
  }
  T val = 0.0;
  for(int i=0;i<val_.size();i++) val += val_[i]*A.val_[i];
  return val;
}


template <typename T>
void tensor<T>::plus(tensor& ret, const tensor& B, T alpha, T beta) const
{
  if( index_ != B.index_ )
  {
    cerr << "Error in tensor plus: dimension doesn't match.\n A: ";
    for(int i=0;i<index_.size();i++)
      cerr << index_[i] << " ";
    cerr << "\n B: ";
    for(int i=0;i<B.index_.size();i++)
      cerr << B.index_[i] << " ";
    cerr << endl;
    return ;
  }
  
  ret.index_ = index_;
  ret.val_.resize(val_.size());
  for(int i=0;i<val_.size();i++)
    ret.val_[i] = alpha * val_[i] + beta * B.val_[i];
}


template <typename T>
void tensor<T>::times(T alpha, T beta) 
{
  for(int i=0;i<val_.size();i++) val_[i] = alpha * val_[i] + beta;
}

  
template <>
void tensor<double>::gemm(double* out, double* in, char transa, char transb,
                          int m, int n, int k, double alpha, double beta)
{
  dgemm(transa, transb, m, n, k, begin(), in, out, alpha, beta);
}

template <>
void tensor<complex<double> >::gemm(complex<double>* out, complex<double>* in,
                               char transa, char transb, int m, int n, int k,
                               complex<double> alpha, complex<double> beta)
{
  zgemm(transa, transb, m, n, k, begin(), in, out, alpha, beta);
}


template <typename T>
void tensor<T>::contract(tensor& A, tensor& B, char transa, char transb,
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
    return ;
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
  A.gemm(begin(), B.begin(), transa, transb, m, n, k, 1, 1);
}


template <typename T>
void tensor<T>::contract(const tensor& A, int a, const tensor& B, int b,
                               bool cover, double sign)
{
  if(A.index_[a]!=B.index_[b])
  {
    cerr << "Error in tensor contract: Dimensions do not match. "
         << A.index_[a] <<"!=" << B.index_[b] << endl;
    return ;
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
        A.val_[i1*mid*ar+k*ar+i2] * B.val_[j1*mid*br+k*br+j2] * sign;
}


template <typename T>
void tensor<T>::contract(T* out, T* in, int row, int col,
                         char transa, char transb, int num,
                         double beta, double alpha)
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
  gemm(out, in, transa, transb, m, n, k, alpha, beta);
}


template <>
void tensor<double>::gesvd(double* in, int row, int col,
                           double* s, double* u, double* v)
{
  dgesvd(in, row, col, s, u, v);
}

template <>
void tensor<complex<double> >::gesvd(complex<double>* in, int row, int col,
                               double* s, complex<double>* u, complex<double>* v)
{
  zgesvd(in, row, col, s, u, v);
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
  gesvd(begin(), left, right, S.begin(), U.begin(), V.begin());
  return ret;
}


template <>
void tensor<double>::heev(int n, double* val, double* vec)
{
  dsyev(begin(), n, val, vec);
}

template <>
void tensor<complex<double> >::heev(int n, double* val, complex<double>* vec)
{
  zheev(begin(), n, val, vec);
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
  heev(d, val, vec);
  return d;
}


/* Explicit instantiation */
template class tensor<double>;
template class tensor< complex<double> >;

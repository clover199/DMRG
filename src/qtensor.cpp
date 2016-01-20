
#include "qtensor.h"
#include "useful.h"

template <typename T>
void qtensor<T>::check_dim_(int n)
{
  for(int i=0;i<index_.size();i++)
  {
    if( dim_[i][sym_[n][i]]==0 )
      dim_[i][sym_[n][i]] = val_[n].dimension(i);
    else if( dim_[i][sym_[n][i]]!=val_[n].dimension(i) )
      cerr << "Error in tensor_quantum check_dim_: "
           << dim_[i][sym_[n][i]] << "!=" << val_[n].dimension(i)
           << endl;
  }
}

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
  
  dim_.resize(index_.size());
  for(int i=0;i<dim_.size();i++) dim_[i].resize(index_[i]);
}


template <typename T>
qtensor<T>::qtensor(T* val, const vector< vector<int> >& dim,
                    const vector< vector<int> >& sym)
{
  index_.resize(dim.size());
  for(int i=0;i<index_.size();i++) index_[i] = dim[i].size();
  dim_ = dim;
  sym_ = sym;
  
  int count = 0;
  val_.resize(sym.size());
  for(int i=0;i<sym.size();i++)
  {
    int add = 1;
    val_[i].index_.resize(index_.size());
    for(int j=0;j<index_.size();j++)
    {
      val_[i].index_[j] = dim[j][sym[i][j]];
      add *= dim[j][sym[i][j]];
    }
    val_[i].val_.resize(add);
    for(int j=0;j<add;j++) val_[i].val_[j] = val[count+j];
    count += add;
  }
}
						
	  
template <typename T>
void qtensor<T>::update(tensor<T>& val,
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

  if( loc.size() == index_.size() )
  {
    bool check = true;
    for(int i=0;i<loc.size();i++) if(loc[i]>=index_[i])
    {
      cerr << "Error in tensor_quantum update: "
              "The input " << i+1 << "th index " << loc[i] 
           << " is out of range 0~" << index_[i]-1 << endl;
      check = false;
    }
    if(check)
    {
      sym_.push_back(loc);
      val_.push_back(val);
      check_dim_(sym_.size()-1);
    }
  }
  else
    cerr << "Error in tensor_quantum update: "
            "The number of index doesn't match. \n";
}


template <typename T>
void qtensor<T>::print() const
{
  if(index_.size()==0) cout << "Empty!" << endl;
  else
  {
    cout << "Dimension for the symmetry sectors: "
         << "index=(" << index_[0];
    for(int i=1;i<index_.size();i++) cout << ", " << index_[i];
    cout << ")\n";
    for(int i=0;i<dim_.size();i++)
    {
      cout << "index=" << i << " : ";
      for(int j=0;j<dim_[i].size();j++) cout << "\t" << dim_[i][j];
      cout << endl;
    }
    cout << "Print all non-zero symmetry sectors: \n";
    for(int i=0;i<val_.size();i++)
    {
      cout << i << " symmetry: (" << sym_[i][0];
      for(int j=1;j<index_.size();j++) cout << ", " << sym_[i][j];
      cout << ")  dimension: (" << val_[i].dimension(0);
      for(int j=1;j<val_[i].index();j++)
        cout << ", " << val_[i].dimension(j);
      cout << ")\n";
    }
  }
}


template <typename T>
void qtensor<T>::print_matrix() const
{
  if(2==index_.size())
  {
    cout << "Print the tensor_quantum as a matrix: \n";
    for(int n=0;n<val_.size();n++)
    {
      cout << "Symmetry: " << sym_[n][0] << " " << sym_[n][1] << endl;
      cout << "\t row=" << val_[n].dimension(0)
           << "  col=" << val_[n].dimension(1) << endl;
      for(int i=0;i<val_[n].dimension(0);i++)
      {
        for(int j=0;j<val_[n].dimension(1);j++)
        {
          T val = val_[n].val_[i*val_[n].dimension(1)+j];
          if(abs(val)>ZERO) cout << val << "\t";
          else cout << (T)0 << "\t";
        }
        cout << endl;
      }
    }
  }
  else cout << "This " << index_.size() << "-tensor cannot be printed "
               "as a matrix\n";
}


template<typename T>
qtensor<T> qtensor<T>::conjugate() const
{
  qtensor<T> ret = *this;
  for(int i=0;i<val_.size();i++) ret.val_[i] = val_[i].conjugate();
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::exchange(int a, int b) const
{
  if( a>=index_.size() or b>=index_.size() or a<0 or b<0)
  {
    cerr << "Error in tensor_quantum exchange: The input index number "
         << a << " and " << b 
         << " is out of range 0~" << index_.size()-1 << endl;
    return *this;
  }
  
  qtensor<T> ret = *this; 
  ret.index_[a] = index_[b];
  ret.index_[b] = index_[a];
  ret.dim_[a] = dim_[b];
  ret.dim_[b] = dim_[a];
  for(int i=0;i<val_.size();i++)
  {
    ret.sym_[i][a] = sym_[i][b];
    ret.sym_[i][b] = sym_[i][a];
    ret.val_[i] = val_[i].exchange(a,b);
  }
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::shift(int num) const
{
  if( num>=index_.size() or num<0)
  {
    cerr << "Error in tensor_quantum transpose: The input index number "
         << num << " is out of range 1~" << index_.size()-1 << endl;
    return *this;
  }
  
  qtensor<T> ret = *this;
  for(int i=0;i<index_.size();i++)
  {
    ret.index_[i] = index_[(i+num)%index_.size()];
    ret.dim_[i] = dim_[(i+num)%index_.size()];
  }
  for(int i=0;i<val_.size();i++)
  {
    for(int j=0;j<index_.size();j++)
      ret.sym_[i] = sym_[(i+num)%index_.size()];
    ret.val_[i] = val_[i].shift(num);
  }
  return ret; 
}


template <typename T>
qtensor<T> qtensor<T>::combine(int min, int max) const
{
  if(min==max) return *this;
  if(min<0 or min>max or max>=index_.size())
  {
    cerr << "Error in tensor_quantum combine: "
            "The imput indexes " << min << " " << max 
         << " is out of range 0~" << index_.size()-1 << endl;
    return *this;
  }
  
  qtensor<T> ret;
  ret.index_.resize(index_.size()-(max-min));
  for(int i=0;i<min;i++) ret.index_[i] = index_[i];
  for(int i=max+1;i<index_.size();i++) ret.index_[i-(max-min)] = index_[i];

  vector< vector<int> > map;
  vector<int> sym_sum;
  sym_sum.resize(max-min+1);
  for(int s=0;s<max-min+1;s++) sym_sum[s] = index_[s+min];
  generate_map(map, sym_sum);
  sym_sum = vector<int> (map.size(),0);
  for(int i=0;i<map.size();i++) for(int s=0;s<max-min+1;s++)
    sym_sum[i] = add(sym_sum[i], map[i][s]);
  int sym = 0;
  for(int i=0;i<sym_sum.size();i++) if(sym<sym_sum[i]) sym = sym_sum[i];
  ret.index_[min] = ++sym;

  ret.dim_.resize(ret.index_.size());
  for(int i=0;i<min;i++) ret.dim_[i] = dim_[i];
  for(int i=max+1;i<index_.size();i++) ret.dim_[i-(max-min)] = dim_[i];
  
  vector<int> sym_dim(map.size(),1);
  for(int i=0;i<map.size();i++) for(int s=0;s<max-min+1;s++)
    sym_dim[i] *= dim_[min+s][map[i][s]];
  ret.dim_[min] = vector<int> (ret.index_[min], 0);
  for(int i=0;i<sym_dim.size();i++) ret.dim_[min][sym_sum[i]] += sym_dim[i];

  generate_map(map, ret.index_);
  vector<bool> check(map.size(), false);
  for(int t=0;t<val_.size();t++)
  {
    vector<int> loc (ret.index_.size(), 0);
    for(int i=0;i<min;i++) loc[i] = sym_[t][i];
    for(int i=min;i<=max;i++) loc[min] = add(loc[min], sym_[t][i]);
    for(int i=max+1;i<index_.size();i++) loc[i-(max-min)] = sym_[t][i];
    int k = 0;
    for(int i=0;i<loc.size();i++) k = k*ret.index_[i] + loc[i];
    check[k] = true;
  }
  int count = 0;  // count the number of different possible symmetry sectors
  for(int i=0;i<check.size();i++) if(check[i]) count++;
  ret.sym_.resize(count);
  ret.val_.resize(count);
  count = 0;  // count the index for sym_ and val_
  // construct and claim space for tensors in val_
  for(int t=0;t<check.size();t++) if(check[t])
  {
    ret.sym_[count] = map[t];
    ret.val_[count].index_.resize(ret.index_.size());
    for(int i=0;i<ret.index_.size();i++)
      ret.val_[count].index_[i] = ret.dim_[i][map[t][i]];
    ret.val_[count].val_.resize(ret.val_[count].dim_());
    count++;
  }
  map.clear();

  // construct the beginning of the symmetry sectors 
  vector<int> store(ret.index_[min],0);
  vector<int> temp(ret.index_[min],0);
  for(int i=0;i<sym_dim.size();i++)
  {
    temp[sym_sum[i]] = sym_dim[i];
    sym_dim[i] = store[sym_sum[i]];
    store[sym_sum[i]] += temp[sym_sum[i]];
  }
  
  for(int t=0;t<val_.size();t++)
  {
    int loc = -1;
    int sum = 0;
    for(int i=min;i<=max;i++) sum = add(sum, sym_[t][i]);
    for(int s=0;s<ret.sym_.size();s++) if(sum == ret.sym_[s][min])
    {
      bool test = true;
      for(int i=0;i<min;i++) if(sym_[t][i] != ret.sym_[s][i]) test = false;
      for(int i=max+1;i<index_.size();i++)
        if(sym_[t][i] != ret.sym_[s][i-(max-min)]) test = false;
      if(test) loc = s;
    }
    if(loc==-1) cerr << "Error in tensor_quantum combine: cannot find sector!\n";
    int beg = 0;
    for(int i=min;i<=max;i++) beg = beg*index_[i] + sym_[t][i];
    beg = sym_dim[beg];

    int l = 1;
    for(int i=0;i<min;i++) l *= val_[t].index_[i];
    int m = 1;
    for(int i=min;i<=max;i++) m *= val_[t].index_[i];
    int r = 1;
    for(int i=max+1;i<index_.size();i++) r *= val_[t].index_[i];
    int size = ret.val_[loc].val_.size();
    for(int i=0;i<l;i++) for(int j=0;j<m;j++) for(int k=0;k<r;k++)
      ret.val_[loc].val_[i*(size/l)+(j+beg)*r+k] = val_[t].val_[i*m*r+j*r+k];
  }
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::split(int index, vector< vector<int> > dim) const
{
  qtensor<T> ret;
  ret.index_.resize(index_.size()+dim.size()-1);
  for(int i=0;i<index;i++) ret.index_[i] = index_[i];
  for(int i=0;i<dim.size();i++) ret.index_[index+i] = dim[i].size();
  for(int i=index+dim.size();i<ret.index_.size();i++)
    ret.index_[i] = index_[i-dim.size()+1];
  
  ret.dim_.resize(index_.size()+dim.size()-1);
  for(int i=0;i<index;i++) ret.dim_[i] = dim_[i];
  for(int i=0;i<dim.size();i++) ret.dim_[index+i] = dim[i];
  for(int i=index+dim.size();i<ret.index_.size();i++)
    ret.dim_[i] = dim_[i-dim.size()+1];

  vector<int> index_new;
  index_new.resize(dim.size());
  for(int i=0;i<index_new.size();i++) index_new[i] = dim[i].size();
  vector< vector<int> > map;
  generate_map(map, index_new);
  vector< vector<int> > sym_dim;
  sym_dim.resize(map.size());
  for(int i=0;i<index_new.size();i++) index_new[i] = 0; //used to store the beg
  for(int i=0;i<map.size();i++)
  {
    sym_dim[i].resize(3);
    sym_dim[i][0] = 0; // the symmetry sum
    for(int s=0;s<dim.size();s++) sym_dim[i][0] = add(sym_dim[i][0], map[i][s]);
    sym_dim[i][1] = 1; // dimension
    for(int s=0;s<dim.size();s++) sym_dim[i][1] *= dim[s][map[i][s]];
    sym_dim[i][2] = index_new[sym_dim[i][0]];
    index_new[sym_dim[i][0]] += sym_dim[i][1];
  }

  index_new.resize(ret.index_.size());
  tensor<T> temp;
  for(int t=0;t<sym_.size();t++) for(int s=0;s<sym_dim.size();s++)
  if(sym_dim[s][0]==sym_[t][index])
  {
    bool check = false;
    int l = 1;
    for(int i=0;i<index;i++) l *= val_[t].index_[i];
    int m = val_[t].index_[index];
    int r = 1;
    for(int i=index+1;i<index_.size();i++) r *= val_[t].index_[i];
    temp.val_.resize(l*sym_dim[s][1]*r);
    int count =0;
    for(int i=0;i<l;i++) for(int j=0;j<sym_dim[s][1];j++) for(int k=0;k<r;k++)
    {
      T store = val_[t].val_[i*m*r+(j+sym_dim[s][2])*r+k];
      if(not check && abs(store)>TOL) check = true;
      temp.val_[count++] = store;
    } 
    if(check)
    {
      for(int i=0;i<index;i++) index_new[i] = sym_[t][i];
      for(int i=0;i<dim.size();i++) index_new[index+i] = map[s][i];
      for(int i=index;i<index_.size();i++) index_new[i+dim.size()] = sym_[t][i];
      ret.sym_.push_back(index_new);
      temp.index_.resize(index_new.size());
      for(int i=0;i<index_new.size();i++)
        temp.index_[i] = ret.dim_[i][index_new[i]];
      ret.val_.push_back(temp);
    }
  }
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::remove_symmetry() const
{
  qtensor<T> ret;
  ret.index_ = vector<int> (index_.size(),1);
  ret.dim_.resize(dim_.size());
  for(int i=0;i<dim_.size();i++) ret.dim_[i] = vector<int> (1,dimension(i));
  ret.sym_.resize(1);
  ret.sym_[0] = vector<int> (index_.size(),0);
  ret.val_.resize(1);
  ret.val_[0].index_.resize(index_.size());
  for(int i=0;i<dim_.size();i++) ret.val_[0].index_[i] = ret.dim_[i][0];
  ret.val_[0].val_ = vector<T> (dimension(),0);
  for(int i=0;i<val_.size();i++)
  {
    vector< vector<int> > map;
    for(int j=0;j<val_[i].dimension();j++)
    {
      generate_map(map, val_[i].index_);
      int loc = 0;
      for(int s=0;s<index_.size();s++)
      {
        int sum = map[j][s];
        for(int k=0;k<sym_[i][s];k++) sum += dim_[s][k];
        loc = loc*ret.dimension(s)+sum;
      }
      ret.val_[0].val_[loc] = val_[i].val_[j];
    }
  }
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::cut(int index, vector<int> cutoff) const
{
  if( index>=index_.size() or index<0)
  {
    cerr << "Error in tensor_quantum cut: The input index number "
         << index << " is out of range 0~" << index_.size()-1 << endl;
    return *this;
  }
  
  if(cutoff.size()!=index_[index])
  {
    cerr << "Error in tensor_quantum cut: The input cut size "
         << cutoff.size() << " doesn't match sector number "
         << index_[index] << endl;
    return *this;
  }
  
  qtensor<T> ret;
  ret.index_ = index_;
  ret.dim_.resize(index_.size());
  for(int i=0;i<index_.size();i++) ret.dim_[i].resize(index_[i]);
  ret.sym_ = sym_;
  ret.val_.resize(val_.size());

  for(int i=0;i<val_.size();i++)
  {
    ret.val_[i] = val_[i].resize(index, cutoff[sym_[i][index]]);
    ret.check_dim_(i);
  }
  return ret;
}


template <typename T>
qtensor<T>& qtensor<T>::plus(const qtensor& A, const qtensor& B,
                              T alpha, T beta)
{
  if( A.index_ != B.index_ )
  {
    cerr << "Error in tensor_quantum plus: "
            "indexes do not match.\n" ;
    return *this;
  }
  if( A.dim_ != B.dim_ )
  {
    cerr << "Error in tensor_quantum plus: "
            "dimensions do not match.\n" ;
    return *this;
  }
  if( A.sym_ != B.sym_ )
  {
    cerr << "Error in tensor_quantum plus: "
            "symmetries do not match.\n" ;
    return *this;
  }
  
  index_ = A.index_;
  dim_ = A.dim_;
  sym_ = A.sym_;
  val_.resize(A.val_.size());
  for(int i=0;i<val_.size();i++) val_[i].plus(A.val_[i], B.val_[i], alpha, beta);
  return *this;
}


template <typename T>
qtensor<T>& qtensor<T>::contract(qtensor<T>& A, qtensor<T>& B,
                                  char transa, char transb, int num)
{
  int indexa = 0;
  if( 'N'==transa or 'n'==transa ) indexa = A.index_.size()-num;
  int indexb = B.index_.size()-num;
  if( 'N'==transb or 'n'==transb ) indexb = 0;
 
  bool check = false;  // check symmetry
  for(int i=0;i<num;i++) if(A.index_[indexa+i]!=B.index_[indexb+i]) check = true;
  if(check)
  {
    cerr << "Error in tensor_quantum contract: The symmetries are not the same.\n"
         << "A: ";
    for(int i=0;i<num;i++) cerr << A.index_[indexa+i] << " ";
    cerr << "\nB: ";
    for(int i=0;i<num;i++) cerr << B.index_[indexb+i] << " ";
    cerr << endl;
    return *this;
  }

  // check dimension
  for(int i=0;i<num;i++) if(A.dim_[indexa+i]!=B.dim_[indexb+i])
  {
    cerr << "Error in tensor_quantum contract: "
            "The dimensions are not the same for i=" << i << ".\n";
    cerr << "A: ";
    for(int j=0;j<A.dim_[indexa+i].size();j++) cerr << A.dim_[indexa+i][j] << " ";
    cerr << "\nB: ";
    for(int j=0;j<B.dim_[indexb+i].size();j++) cerr << B.dim_[indexb+i][j] << " ";
    cerr << endl;
    return *this;
  }
  
  int add = 0;
  index_.resize(A.index_.size()+B.index_.size()-2*num);
  for(int i=0;i<A.index_.size();i++) if(i<indexa or i>=indexa+num)
    index_[add++] = A.index_[i];
  for(int i=0;i<B.index_.size();i++) if(i<indexb or i>=indexb+num)
    index_[add++] = B.index_[i];

  dim_.resize(index_.size());
  for(int i=0;i<index_.size();i++) dim_[i] = vector<int>(index_[i],0);

  int ab = A.val_.size()*B.val_.size();
  sym_.resize(ab);
  val_.clear();
  val_.resize(ab);
  vector<int> map(ab,0);  // used to check if the symmetry sector exists.
  int count = 0;
  for(int i=0;i<A.val_.size();i++) for(int j=0;j<B.val_.size();j++)
  {
    check = true;
    for(int k=0;k<num;k++) if(A.sym_[i][indexa+k]!= B.sym_[j][indexb+k])
      check = false;

    if(check)
    {
      add = 0;
      sym_[count].resize(index_.size());
      for(int k=0;k<A.index_.size();k++) if(k<indexa or k>=indexa+num)
        sym_[count][add++] = A.sym_[i][k];
      for(int k=0;k<B.index_.size();k++) if(k<indexb or k>=indexb+num)
        sym_[count][add++] = B.sym_[j][k];
      for(int s=0;s<index_.size();s++)
        map[count] = map[count]*index_[s]+sym_[count][s];
      bool same = false;
      for(int t=0;t<count;t++) if(map[count]==map[t])
      {
        val_[t].contract(A.val_[i], B.val_[j], transa, transb, num);
        same = true;
      }
      if(not same)
      {
        val_[count].contract(A.val_[i], B.val_[j], transa, transb, num);
        check_dim_(count++);
      }
    }
  }
  sym_.resize(count);
  val_.resize(count);
  return *this;
}


template <typename T>
qtensor<T>& qtensor<T>::contract(const qtensor<T>& A, int a,
                                 const qtensor<T>& B, int b)
{
  if(A.index_[a]!=B.index_[b])
  {
    cerr << "Error in tensor_quantum contract: symmetry sector doesn't match. "
         << A.index_[a] << "!=" << B.index_[b] << endl;
    return *this;
  }
  for(int i=0;i<A.index_[a];i++) if(A.dim_[a][i]!=B.dim_[b][i])
  {
    cerr << "Error in tensor_quantum contract: "
            "dimension of symmetry sector " << i << " doesn't match. "
         << A.dim_[a][i] << "!=" << B.dim_[b][i] << endl;
    return *this;
  }
  
  index_.resize(A.index_.size()+B.index_.size()-2);
  for(int i=0;i<a;i++) index_[i] = A.index_[i];
  for(int i=a;i<a+b;i++) index_[i] = B.index_[i-a];
  for(int i=a+b;i<a+B.index_.size()-1;i++)
    index_[i] = B.index_[i-a+1];
  for(int i=a+B.index_.size()-1;i<index_.size();i++)
    index_[i] = A.index_[i-B.index_.size()+2];

  dim_.resize(index_.size());
  for(int i=0;i<index_.size();i++) dim_[i] = vector<int>(index_[i],0);

  int ab = A.val_.size()*B.val_.size();
  sym_.resize(ab);
  val_.clear();
  val_.resize(ab);
  vector<int> map(ab,0);
  int count=0;
  for(int l=0;l<A.val_.size();l++) for(int r=0;r<B.val_.size();r++)
  if(A.sym_[l][a]==B.sym_[r][b])
  {
    sym_[count].resize(index_.size());
    for(int i=0;i<a;i++) sym_[count][i] = A.sym_[l][i];
    for(int i=a;i<a+b;i++) sym_[count][i] = B.sym_[r][i-a];
    for(int i=a+b;i<a+B.index_.size()-1;i++)
      sym_[count][i] = B.sym_[r][i-a+1];
    for(int i=a+B.index_.size()-1;i<index_.size();i++)
      sym_[count][i] = A.sym_[l][i-B.index_.size()+2];
    for(int s=0;s<index_.size();s++)
      map[count] = map[count]*index_[s]+sym_[count][s];
    bool same = false;
    for(int t=0;t<count;t++) if(map[count]==map[t])
    {
      val_[t].contract(A.val_[l], a, B.val_[r], b);
      same = true;
    }
    if(not same)
    {
      val_[count].contract(A.val_[l], a, B.val_[r], b);
      check_dim_(count++);
    }
  }
  sym_.resize(count);
  val_.resize(count);
  for(int i=0;i<val_.size();i++) check_dim_(i);
  return *this;
}


template <typename T>
void qtensor<T>::contract(T* out, T* in,
                          const vector< vector<int> >& map,
                          char transa, char transb, int num)
{
  for(int i=0;i<map.size();i++)
    val_[map[i][0]].contract(out+map[i][2], in+map[i][1],
                     map[i][3], map[i][4], transa, transb, num, map[i][5]);
}


template <typename T>
bool qtensor<T>::get_map(vector< vector<int> >& map_ret,
                         vector< vector<int> >& dim_ret,
                         vector< vector<int> >& sym_ret,
                         const vector< vector<int> >& dim,
                         const vector< vector<int> >& sym,
                         char transa, char transb, int num, double beta) const
{
  int indexa = 0;
  if( 'N'==transa or 'n'==transa ) indexa = index_.size()-num;
  int indexb = dim.size()-num;
  if( 'N'==transb or 'n'==transb ) indexb = 0;

  bool check = false;  // check symmetry
  for(int i=0;i<num;i++) if(index_[indexa+i]!=dim[indexb+i].size()) check = true;
  if(check)
  {
    cerr << "Error in tensor_quantum get_map: The symmetries are not the same.\n"
         << "Tensor: ";
    for(int i=0;i<num;i++) cerr << index_[indexa+i] << " ";
    cerr << "\nArray:  ";
    for(int i=0;i<num;i++) cerr << dim[indexb+i].size() << " ";
    cerr << endl;
    return false;
  }

  // check dimension
  for(int i=0;i<num;i++) if(dim_[indexa+i]!=dim[indexb+i])
  {
    cerr << "Error in tensor_quantum get_map: "
            "The dimensions are not the same for i=" << i << ".\n";
    cerr << "tensor: ";
    for(int j=0;j<dim_[indexa+i].size();j++) cerr << dim_[indexa+i][j] << " ";
    cerr << "\narray:  ";
    for(int j=0;j<dim[indexb+i].size();j++) cerr << dim[indexb+i][j] << " ";
    cerr << endl;
    return false;
  }

  int add = 0;  // get index
  vector<int> index;
  index.resize(index_.size()+dim.size()-2*num, 0);
  for(int i=0;i<index_.size();i++) if(i<indexa or i>=indexa+num)
    index[add++] = index_[i];
  for(int i=0;i<dim.size();i++) if(i<indexb or i>=indexb+num)
    index[add++] = dim[i].size();

  add = 0;  // get dim_ret
  dim_ret.resize(index.size());
  for(int i=0;i<index_.size();i++) if(i<indexa or i>=indexa+num)
    dim_ret[add++] = dim_[i];
  for(int i=0;i<dim.size();i++) if(i<indexb or i>=indexb+num)
    dim_ret[add++] = dim[i];

  vector< vector<int> > sym_dim;  // dim, dim of row, begin point of sym
  sym_dim.resize(sym.size());
  for(int i=0;i<sym.size();i++)
  {
    sym_dim[i] = vector<int> (3,1);
    for(int s=0;s<dim.size();s++) sym_dim[i][0] *= dim[s][sym[i][s]];
    if(indexb) for(int s=0;s<indexb;s++) sym_dim[i][1] *= dim[s][sym[i][s]];
    else for(int s=0;s<num;s++) sym_dim[i][1] *= dim[s][sym[i][s]];
  }
  sym_dim[0][2] = 0;
  for(int i=1;i<sym.size();i++) sym_dim[i][2] = sym_dim[i-1][2] + sym_dim[i-1][0];

  int ab = val_.size()*sym.size();
  map_ret.resize(ab);
  sym_ret.resize(ab);
  vector<int> map(ab,0);  // used to check if the symmetry sector exists.
  vector<int> map_dim(ab+1,0); // used to store the beginning of the sym_ret
  int count_map = 0;
  int count_sym = 0;
  for(int i=0;i<val_.size();i++) for(int j=0;j<sym.size();j++)
  {
    check = true;
    for(int k=0;k<num;k++) if(sym_[i][indexa+k]!= sym[j][indexb+k])
      check = false;
    if(check)
    {
      sym_ret[count_sym].resize(index.size());
      add = 0;
      for(int k=0;k<index_.size();k++) if(k<indexa or k>=indexa+num)
        sym_ret[count_sym][add++] = sym_[i][k];
      for(int k=0;k<dim.size();k++) if(k<indexb or k>=indexb+num)
        sym_ret[count_sym][add++] = sym[j][k];
      for(int s=0;s<index.size();s++)
        map[count_sym] = map[count_sym]*index[s]+sym_ret[count_sym][s];
      map_ret[count_map].resize(6);
      map_ret[count_map][0] = i;
      map_ret[count_map][1] = sym_dim[j][2];
      map_ret[count_map][2] = map_dim[count_sym];
      map_ret[count_map][3] = sym_dim[j][1];
      map_ret[count_map][4] = sym_dim[j][0]/sym_dim[j][1];
      map_ret[count_map][5] = beta;
      bool same = false;
      for(int t=0;t<count_sym;t++) if(map[count_sym]==map[t])
      {
        map_ret[count_map][2] = map_dim[t];
        map_ret[count_map][5] = 1;
        same = true;
      }
      if(not same)
      {
        map_dim[count_sym+1] = 1;
        for(int s=0;s<index.size();s++)
          map_dim[count_sym+1] *= dim_ret[s][sym_ret[count_sym][s]];
        map_dim[count_sym+1] += map_dim[count_sym];
        count_sym++;
      }
      count_map++;
    }
  }
  map_ret.resize(count_map);
  sym_ret.resize(count_sym);
  return true;  
}


template <typename T>
vector<double> qtensor<T>::svd(qtensor& U, qtensor<double>& S, qtensor& V,
                               int num, int cutoff)
{
  if( num>=index_.size() or num<=0 )
  {
    cerr << "Error in tensor_quantum svd: "
            "The number of left indexes should be 1~"<< index_.size()-1 << endl;
    return vector<double>();
  }
  
  qtensor<T> comb;
  comb = combine(num, index_.size()-1);
  comb = comb.combine(0, num-1);

  int sym = add(comb.sym_[0][0], comb.sym_[0][1]);
  bool check = true;
  for(int i=0;i<comb.sym_.size();i++)
    if(add(comb.sym_[i][0], comb.sym_[i][1])!=sym) check = false;
  if(not check)
  {
    cerr << "Warning in tensor_quantum svd: "
            "symmetry sectors do not have constraint.\n";
    comb = comb.remove_symmetry();
    return comb.svd(U, S, V, num, cutoff);
  }
  
  U.index_ = comb.index_;
  U.index_[1] = U.index_[0];
  V.index_ = comb.index_;
  V.index_[0] = V.index_[1];
  S.index_ = comb.index_;

  U.dim_ = comb.dim_;
  U.dim_[1] = U.dim_[0];
  V.dim_ = comb.dim_;
  V.dim_[0] = V.dim_[1];
  S.dim_ = comb.dim_;

  sym = comb.val_.size();
  U.sym_.resize(sym);
  V.sym_.resize(sym);
  S.sym_.resize(sym);
  U.val_.resize(sym);
  V.val_.resize(sym);
  S.val_.resize(sym);
  int sum = 0; // number of non-zero singular values
  for(int i=0;i<sym;i++)  // svd for different symmetry sectors
  {
    comb.val_[i].svd(U.val_[i], S.val_[i], V.val_[i]);
    sum += S.val_[i].val_.size();
    U.sym_[i] = comb.sym_[i];
    U.sym_[i][1] = U.sym_[i][0];
    V.sym_[i] = comb.sym_[i];
    V.sym_[i][0] = V.sym_[i][1];
    S.sym_[i] = comb.sym_[i];
  }

  vector<double> s; // store the sorted singular values
  s.resize(sum);
  // set the cutoff for different symmetry sectors
  sum = sum < cutoff ? sum : cutoff;
  vector<int> u_num(sym,0);  // the cutoff for different sectors of U
  for(int t=0;t<sum;t++)
  {
    double max = -1;
    int loc;
    for(int i=0;i<sym;i++)
    if(u_num[i]<S.val_[i].val_.size() && max<S.val_[i].val_[u_num[i]])
    {
      max = S.val_[i].val_[u_num[i]];
      loc = i;
    }
    s[t] = max;
    u_num[loc] += 1;
  }
  vector<int> v_num(u_num);  // the cutoff for different sectors for V
  for(int t=sum;t<s.size();t++)
  {
    double max = -1;
    int loc;
    for(int i=0;i<sym;i++)
    if(u_num[i]<S.val_[i].val_.size() && max<S.val_[i].val_[u_num[i]])
    {
      max = S.val_[i].val_[u_num[i]];
      loc = i;
    }
    s[t] = max;
    u_num[loc] += 1;
  }
  u_num = v_num;

  for(int i=0;i<sym;i++)
  {
    S.val_[i].diag(2);
    if(comb.val_[i].index_[0]<comb.val_[i].index_[1])
      S.val_[i] = S.val_[i].resize(1,comb.val_[i].index_[1]);
    else S.val_[i] = S.val_[i].resize(0,comb.val_[i].index_[0]);
  }
  if(cutoff==0) return s;

  // if cutoff is larger than the number of non-zero singular values
  // adding the zeros for U
  int diff = comb.dimension(0);  // dimension
  diff = diff < cutoff ? diff : cutoff;  // the actuall cutoff
  diff = diff-sum; 
  while(diff>0) for(int i=0;i<sym;i++) if(u_num[i]<comb.dim_[0][i] && diff>0)
  {
    u_num[i] += 1;
    diff--;
  }
  // adding the zeros for V
  diff = comb.dimension(1);
  diff = diff < cutoff ? diff : cutoff;
  diff = diff-sum;
  while(diff>0) for(int i=0;i<sym;i++) if(v_num[i]<comb.dim_[1][i] && diff>0)
  {
    v_num[i] += 1;
    diff--;
  }
  
  for(int i=0;i<sym;i++)
  {
    U.val_[i] = U.val_[i].resize(1, u_num[i]);
    U.dim_[1][U.sym_[i][1]] = u_num[i];
    V.val_[i] = V.val_[i].resize(0, v_num[i]);
    V.dim_[0][V.sym_[i][0]] = v_num[i];
    S.val_[i] = S.val_[i].resize(0, u_num[i]);
    S.dim_[0][S.sym_[i][0]] = u_num[i];
    S.val_[i] = S.val_[i].resize(1, v_num[i]);
    S.dim_[1][S.sym_[i][1]] = v_num[i];
  }

  vector< vector<int> > dim;
  dim.resize(num);
  for(int i=0;i<num;i++) dim[i] = dim_[i];
  U = U.split(0, dim);
  dim.resize(index_.size()-num);
  for(int i=0;i<dim.size();i++) dim[i] = dim_[i+num];
  V = V.split(1, dim);
  return s;
}


template <typename T>
int qtensor<T>::eig(double* val, T* vec, int sector)
{
  if(index_.size()!=2)
  {
    cerr << "Error in tensor_quantum eig: cannot diagonalize. "
            "This is a " << index_.size() << "-tensor.";
    return 0;
  }
  if(index_[0]!=index_[1])
  {
    cerr << "Error in tensor_quantum eig: the symmetries are not the same.\n"
         << index_[0] << "!=" << index_[1] << endl;
    return 0;
  }
  for(int i=0;i<index_[0];i++) if(dim_[0][i]!=dim_[1][i])
  {
    cerr << "Error in tensor_quantum eig: the symmetry sector " << i
         << " is not square.\n"
         << dim_[i][0] << "!=" << dim_[i][1] << endl;
    return 0;
  }
  bool check = false;
  for(int i=0;i<val_.size();i++) if(sym_[i][0]!=sym_[i][1])
  {
    cerr << "Warning in tensor_quantum eig: not block diagonalized."
            "Combining the symmetry sectors\n";
    check = true;
    return 0;
  }
  if(sym_.size()==1) return val_[0].eig(val, vec);
  if(check or sector==-1)
  {
    qtensor<T> temp;
    temp = remove_symmetry();
    return temp.eig(val, vec, -1);
  }
  int loc = -1;
  for(int i=0;i<sym_.size();i++) if(sym_[i][0]==sector) loc = i;
  if(loc==-1)
  {
    cerr << "Error in tensor_quantum eig: cannot find symmetry sector "
         << sector << endl;
    return 0;
  }
  return val_[loc].eig(val, vec);
}


template class qtensor<double>;
template class qtensor< complex<double> >;

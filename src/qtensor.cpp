
#include "qtensor.h"
#include "useful.h"

template <typename T>
qtensor<T>::qtensor(int i0, int i1, int i2, int i3,
                    int i4, int i5, int i6, int i7)
{
  vector<int> index; index.clear();
  if(i0>0){ index.push_back(i0);
  if(i1>0){ index.push_back(i1);
  if(i2>0){ index.push_back(i2);
  if(i3>0){ index.push_back(i3);
  if(i4>0){ index.push_back(i4);
  if(i5>0){ index.push_back(i5);
  if(i6>0){ index.push_back(i6);
  if(i7>0){ index.push_back(i7); }}}}}}}}

  dim_.resize(index.size());
  for(int i=0;i<dim_.size();i++) dim_[i].resize(index[i]);
  dir_.resize(dim_.size());
}


template <typename T>
qtensor<T>::qtensor(T* val, const vector< vector<int> >& dim,
                    const vector< vector<int> >& sym)
{
  dir_ = vector<int> (dim.size(),0);
  dim_ = dim;
  sym_ = sym;

  int count = 0;
  val_.resize(sym.size());
  for(int i=0;i<sym.size();i++)
  {
    int d = 1;
    val_[i].index_.resize(dim.size());
    for(int j=0;j<dim.size();j++)
    {
      val_[i].index_[j] = dim[j][sym[i][j]];
      d *= dim[j][sym[i][j]];
    }
    val_[i].val_.resize(d);
    for(int j=0;j<d;j++) val_[i].val_[j] = val[count+j];
    count += d;
  }
}


template <typename T>
void qtensor<T>::update(tensor<T>& val,
                        int i0, int i1, int i2, int i3,
                        int i4, int i5, int i6, int i7)
{
  if(val.empty()) return;

  vector<int> loc;
  if(i0>=0){ loc.push_back(i0);
  if(i1>=0){ loc.push_back(i1);
  if(i2>=0){ loc.push_back(i2);
  if(i3>=0){ loc.push_back(i3);
  if(i4>=0){ loc.push_back(i4);
  if(i5>=0){ loc.push_back(i5);
  if(i6>=0){ loc.push_back(i6);
  if(i7>=0){ loc.push_back(i7);}}}}}}}}

  if( loc.size() != dim_.size() )
  {
    cerr << "Error in tensor_quantum update: "
            "The number of index doesn't match. \n";
    return ;
  }

  for(int i=0;i<loc.size();i++) if(loc[i]>=dim_[i].size())
  {
    cerr << "Error in tensor_quantum update: "
            "The input " << i+1 << "th index " << loc[i]
         << " is out of range 0~" << dim_[i].size()-1 << endl;
    return ;
  }

  sym_.push_back(loc);
  val_.push_back(val);
  check_dim_(sym_.size()-1);
}


template <typename T>
void qtensor<T>::u_dim(int i0, int i1, int i2, int i3,
                       int i4, int i5, int i6, int i7)
{
  vector<int> index; index.clear();
  if(i0>0){ index.push_back(i0);
  if(i1>0){ index.push_back(i1);
  if(i2>0){ index.push_back(i2);
  if(i3>0){ index.push_back(i3);
  if(i4>0){ index.push_back(i4);
  if(i5>0){ index.push_back(i5);
  if(i6>0){ index.push_back(i6);
  if(i7>0){ index.push_back(i7); }}}}}}}}

  if( index.size() != dim_.size() )
  {
    cerr << "Error in tensor_quantum u_dim: "
            "The number of index doesn't match. \n";
    return;
  }

  for(int i=0;i<index.size();i++) for(int j=0;j<dim_[i].size();j++)
  {
    if( dim_[i][j]!=0 and dim_[i][j]!=index[i])
    {
      cerr << "Error in tensor_quantum u_dim: "
              "The input " << i+1 << "th index " << index[i]
           << " doesn't match existing " << dim_[i][j]
           << " of symmetry sector " << j << endl;
      return;
    }
  }

  for(int i=0;i<index.size();i++) for(int j=0;j<dim_[i].size();j++)
    dim_[i][j] = index[i];
}


template <typename T>
qtensor<T> qtensor<T>::id() const
{
  qtensor<T> ret;
  if(dim_[0]!=dim_[2] or dim_[1]!=dim_[3])
  {
    cerr << "Error in tensor_quantum id: Wrong dim" << endl;
    return ret;
  }
  ret.dim_ = dim_;
  ret.dir_ = dir_;
  ret.sym_.resize(dim_[0].size()*dim_[1].size());
  ret.val_.resize(dim_[0].size()*dim_[1].size());

  for(int s1=0;s1<dim_[0].size();s1++) for(int s2=0;s2<dim_[1].size();s2++)
  {
    int index = s1*dim_[1].size()+s2;
    vector<int> sym(4,0);
	sym[0]=s1; sym[1]=s2; sym[2]=s1; sym[3]=s2;
    ret.sym_[index] = sym;
	vector<int> dim;
	get_index_(dim, sym);
    ret.val_[index] = tensor<T> (dim);
	ret.val_[index].id(0,2, 1,3);
  }
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::id(const qtensor<T>& bulk_mpo) const
{
  qtensor<T> ret;
  ret.dim_ = bulk_mpo.dim_;
  ret.dim_[0] = dim_[1];
  ret.dim_[2] = dim_[1];

  ret.dir_ = bulk_mpo.dir_;
  ret.dir_[0] = dir_[1];
  ret.dir_[2] = dir_[1];

  ret.sym_.clear();
  ret.val_.clear();
  vector<int> sym(ret.dim_.size(), 0);
  tensor<T> temp;
  vector< vector<int> > map;
  for(int s0=0;s0<ret.dim_[0].size();s0++) for(int s1=0;s1<ret.dim_[1].size();s1++)
  {
    sym[0] = s0; sym[1] = s1;
    sym[2] = s0; sym[3] = s1;
    ret.sym_.push_back(sym);

    temp = tensor<T> (ret.dim_[0][s0], ret.dim_[1][s1], ret.dim_[0][s0], ret.dim_[1][s1]);
    generate_map(map, temp.index_);
    for(int i=0;i<map.size();i++) if(map[i][0]==map[i][2] and map[i][1]==map[i][3])
      temp.val_[i] = 1;
    ret.val_.push_back(temp);
  }
  return ret;
}


template <typename T>
void qtensor<T>::add_sign(int i0, int i1, int i2, int i3,
                          int i4, int i5, int i6, int i7)
{
  dir_.clear();
  if(i0!=2){ dir_.push_back(i0);
  if(i1!=2){ dir_.push_back(i1);
  if(i2!=2){ dir_.push_back(i2);
  if(i3!=2){ dir_.push_back(i3);
  if(i4!=2){ dir_.push_back(i4);
  if(i5!=2){ dir_.push_back(i5);
  if(i6!=2){ dir_.push_back(i6);
  if(i7!=2){ dir_.push_back(i7);}}}}}}}}

  if( dir_.size() != dim_.size() )
    cerr << "Error in tensor_quantum add_sign: "
            "The number of index doesn't match. \n";
}


template <typename T>
void qtensor<T>::print() const
{
  if(dim_.size()==0) cout << "Empty!" << endl;
  else
  {
    cout << "This is a fermionic tensor with sigh direction: (" << dir_[0];
    for(int i=1;i<dir_.size();i++) cout << ", " << dir_[i];
    cout << ")\n";

    cout << "Dimension for different symmetry sectors:\n";
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
      for(int j=1;j<dim_.size();j++) cout << ", " << sym_[i][j];
      cout << ")  dimension: (" << val_[i].dimension(0);
      for(int j=1;j<val_[i].index();j++)
        cout << ", " << val_[i].dimension(j);
      cout << ")\n";
    }
  }
}


template <typename T>
void qtensor<T>::print_all(bool all) const
{
  print();
  if(val_.size()!=0) for(int n=0;n<val_.size();n++)
  {
    cout << "Print all non-zero values of the tensor " << n << ":\n";
    vector< vector<int> > map;
    generate_map(map, val_[n].index_);
    for(int i=0;i<val_[n].val_.size();i++) if(all)
    {
      cout << "(" << map[i][0];
      for(int j=1;j<val_[n].index_.size();j++) cout << ", " << map[i][j];
      cout << ") " << val_[n].val_[i] << "\n";
    }
    else if(abs(val_[n].val_[i])>ZERO)
    {
      cout << "(" << map[i][0];
      for(int j=1;j<val_[n].index_.size();j++) cout << ", " << map[i][j];
      cout << ") " << val_[n].val_[i] << "\n";
    }
  }
}


template <typename T>
void qtensor<T>::print_matrix() const
{
  cout.precision(3);
  if(2==dim_.size())
  {
    cout << "This is a fermion operator. dir_=("
         << dir_[0] << ", " << dir_[1] << ")\n";

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
  else cout << "This " << dim_.size() << "-tensor cannot be printed "
               "as a matrix\n";
}


template<typename T>
qtensor<T> qtensor<T>::conjugate()
{
  for(int i=0;i<dir_.size();i++) dir_[i] *= -1;
  for(int i=0;i<val_.size();i++) val_[i].conjugate();
  return *this;
}


template <>
qtensor< complex<double> > qtensor<double>::comp() const
{
  qtensor<complex<double> > ret;
  ret.dim_ = dim_;
  ret.sym_ = sym_;
  ret.dir_ = dir_;

  ret.val_.resize(val_.size());
  for(int i=0;i<val_.size();i++) val_[i].to_complex(ret.val_[i]);
  return ret;
}


template <>
qtensor< complex<double> > qtensor< complex<double> >::comp() const
{
  return *this;
}

template <typename T>
qtensor<T> qtensor<T>::simplify() const
{
  qtensor<T> ret;
  ret.dir_ = dir_;
  ret.dim_ = dim_;
  ret.sym_.clear();
  ret.val_.clear();
  for(int i=0;i<sym_.size();i++)
  {
    bool check = false;
    for(int j=0;j<val_[i].val_.size();j++)
      if(abs(val_[i].val_[j])>TOL) check = true;
    if(check)
    {
      ret.sym_.push_back(sym_[i]);
      ret.val_.push_back(val_[i]);
    }
  }
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::exchange(int a, int b) const
{
  if(dim_.size()==0) return *this;
  if( a>=dim_.size() or b>=dim_.size() or a<0 or b<0)
  {
    cerr << "Error in tensor_quantum exchange: The input index number "
         << a << " and " << b
         << " is out of range 0~" << dim_.size()-1 << endl;
    return *this;
  }

  qtensor<T> ret = *this;
  ret.dim_[a] = dim_[b];
  ret.dim_[b] = dim_[a];
  ret.dir_[a] = dir_[b];
  ret.dir_[b] = dir_[a];
  for(int i=0;i<val_.size();i++)
  {
    ret.sym_[i][a] = sym_[i][b];
    ret.sym_[i][b] = sym_[i][a];
    val_[i].exchange(ret.val_[i], a,b);
  }
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::shift(int num) const
{
  if( num>=dim_.size() or num<0)
  {
    cerr << "Error in tensor_quantum transpose: The input index number "
         << num << " is out of range 1~" << dim_.size()-1 << endl;
    return *this;
  }

  qtensor<T> ret = *this;
  for(int i=0;i<dim_.size();i++)
    ret.dim_[i] = dim_[(i+num)%dim_.size()];
  for(int i=0;i<dir_.size();i++)
    ret.dir_[i] = dir_[(i+num)%dir_.size()];
  for(int i=0;i<val_.size();i++)
  {
    for(int j=0;j<dim_.size();j++)
      ret.sym_[i] = sym_[(i+num)%dim_.size()];
    val_[i].shift(ret.val_[i], num);
  }
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::take(int n, int k) const
{
  qtensor<T> ret;
  if(n<0 or n>=dim_.size())
  {
    cout << "Error in tensor_quantum take: "
            "the index " << n << " is out of the range 0~"
         << dim_.size()-1 << endl;
    return ret;
  }

  if(dim_[n].size()!=1)
  {
    cout << "Error in tensor_quantum take: "
            "this index has symmetry " << dim_[n].size() << endl;
    return ret;
  }

  if( k<0 or k>=dim_[n][0])
  {
    cout << "Error in tensor_quantum take: "
            "the number " << k << " is out of the dimension range 0~"
         << dim_[n][0] << endl;
    return ret;
  }

  ret.dim_.resize(dim_.size()-1);
  for(int i=0;i<n;i++) ret.dim_[i] = dim_[i];
  for(int i=n;i<ret.dim_.size();i++) ret.dim_[i] = dim_[i+1];

  ret.dir_.resize(dir_.size()-1);
  for(int i=0;i<n;i++) ret.dir_[i] = dir_[i];
  for(int i=n;i<ret.dir_.size();i++) ret.dir_[i] = dir_[i+1];

  ret.sym_.resize(sym_.size());
  ret.val_.resize(val_.size());
  vector< vector<int> > map;
  for(int s=0;s<sym_.size();s++)
  {
    ret.sym_[s].resize(ret.dim_.size());
    for(int i=0;i<n;i++) ret.sym_[s][i] = sym_[s][i];
    for(int i=n;i<ret.dim_.size();i++) ret.sym_[s][i] = sym_[s][i+1];
    val_[s].take(ret.val_[s], n,k);
  }
  return ret.simplify();
}

template <typename T>
qtensor<T> qtensor<T>::left() const
{
  qtensor<T> ret = take(0, dim_[0][0]-1);
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::right() const
{
  qtensor<T> ret = take(2, 0);
  ret = ret.exchange(0,1);
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::combine(int min, int max) const
{
  if(min==max) return *this;
  if(min<0 or min>max or max>=dim_.size())
  {
    cerr << "Error in tensor_quantum combine: "
            "The input indexes " << min << " " << max
         << " is out of range 0~" << dim_.size()-1 << endl;
    return *this;
  }

  qtensor<T> ret;
  for(int i=min;i<max;i++) if(dir_[i]!=dir_[max])
  {
    cerr << "Error in tensor_quantum combine: "
            "The combining indexes are in different directions: ";
    for(int j=min;j<=max;j++) cerr << dir_[j] << " ";
    cerr << endl;
	return *this;
  }
  ret.dir_.resize(dir_.size()-(max-min));
  for(int i=0;i<min;i++) ret.dir_[i] = dir_[i];
  for(int i=max;i<dir_.size();i++) ret.dir_[i-(max-min)] = dir_[i];

  vector< vector<int> > map;  // all the possible symmetry combinations
  vector<int> sym_sum;
  sym_sum.resize(max-min+1);  // used as index for region min~max
  for(int s=0;s<max-min+1;s++) sym_sum[s] = dim_[s+min].size();
  generate_map(map, sym_sum);
  // calculate the dimension for each symmetry possibilities in map
  vector<int> sym_dim(map.size(),1);
  for(int i=0;i<map.size();i++) for(int s=0;s<max-min+1;s++)
    sym_dim[i] *= dim_[min+s][map[i][s]];
  // calculate the total symmetry for all possibilities in map
  sym_sum = vector<int> (map.size(),0);
  for(int i=0;i<map.size();i++) for(int s=0;s<max-min+1;s++)
    sym_sum[i] = ss(sym_sum[i], map[i][s]);
  int sym = 0;  // find symmetry number for the combined part
  for(int i=0;i<sym_sum.size();i++) if(sym_sum[i]>sym) sym = sym_sum[i];


  ret.dim_.resize(dim_.size()-(max-min));  // get dim_ for ret
  for(int i=0;i<min;i++) ret.dim_[i] = dim_[i];
  for(int i=max+1;i<dim_.size();i++) ret.dim_[i-(max-min)] = dim_[i];
  ret.dim_[min] = vector<int> (sym+1, 0);
  for(int i=0;i<sym_dim.size();i++) ret.dim_[min][sym_sum[i]] += sym_dim[i];

  vector<int> index(ret.dim_.size(), 0);  // index for ret
  for(int i=0;i<index.size();i++) index[i] = ret.dim_[i].size();
  generate_map(map, index);  // all possible symmetries for ret
  vector<bool> check(map.size(), false);
  for(int t=0;t<val_.size();t++)
  {
    vector<int> loc (index.size(), 0);
    for(int i=0;i<min;i++) loc[i] = sym_[t][i];
    for(int i=min;i<=max;i++) loc[min] = ss(loc[min], sym_[t][i]);
    for(int i=max+1;i<dim_.size();i++) loc[i-(max-min)] = sym_[t][i];
    int k = 0;
    for(int i=0;i<loc.size();i++) k = k*index[i] + loc[i];
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
    ret.val_[count].index_.resize(index.size());
    for(int i=0;i<index.size();i++)
      ret.val_[count].index_[i] = ret.dim_[i][map[t][i]];
    ret.val_[count].val_.resize(ret.val_[count].dim_());
    count++;
  }
  map.clear();

  // construct the beginning of the symmetry sectors in sym_dim
  vector<int> store(index[min],0);  // store the beginning for each sector
  for(int i=0;i<sym_dim.size();i++)
  {
    int temp = sym_dim[i];
    sym_dim[i] = store[sym_sum[i]];
    store[sym_sum[i]] += temp;
  }

  for(int t=0;t<val_.size();t++)
  {
    int loc = -1;  // location in ret.val_
    int sum = 0;  // sum of symmetry
    for(int i=min;i<=max;i++) sum = ss(sum, sym_[t][i]);
    for(int s=0;s<ret.sym_.size();s++) if(sum == ret.sym_[s][min])
    {
      bool test = true;
      for(int i=0;i<min;i++)
        if(sym_[t][i] != ret.sym_[s][i]) test = false;
      for(int i=max+1;i<dim_.size();i++)
        if(sym_[t][i] != ret.sym_[s][i-(max-min)]) test = false;
      if(test) loc = s;
    }
    if(loc==-1)
      cerr << "Error in tensor_quantum combine: cannot find sector!\n";
    int beg = 0;
    for(int i=min;i<=max;i++) beg = beg*dim_[i].size() + sym_[t][i];
    beg = sym_dim[beg];

    int l = 1;
    for(int i=0;i<min;i++) l *= val_[t].index_[i];
    int m = 1;
    for(int i=min;i<=max;i++) m *= val_[t].index_[i];
    int r = 1;
    for(int i=max+1;i<dim_.size();i++) r *= val_[t].index_[i];
    int size = ret.val_[loc].val_.size();
    for(int i=0;i<l;i++) for(int j=0;j<m;j++) for(int k=0;k<r;k++)
      ret.val_[loc].val_[i*(size/l)+(j+beg)*r+k] = val_[t].val_[i*m*r+j*r+k];
  }
  return ret;
}


template <typename T>
qtensor<T> qtensor<T>::split(int n, vector< vector<int> > dim) const
{
  qtensor<T> ret;

  ret.dir_.resize(dir_.size()+dim.size()-1);
  for(int i=0;i<n;i++) ret.dir_[i] = dir_[i];
  for(int i=0;i<dim.size();i++) ret.dir_[n+i] = dir_[n];
  for(int i=n+dim.size();i<ret.dir_.size();i++)
    ret.dir_[i] = dir_[i-dim.size()+1];

  ret.dim_.resize(dim_.size()+dim.size()-1);
  for(int i=0;i<n;i++) ret.dim_[i] = dim_[i];
  for(int i=0;i<dim.size();i++) ret.dim_[n+i] = dim[i];
  for(int i=n+dim.size();i<ret.dim_.size();i++)
    ret.dim_[i] = dim_[i-dim.size()+1];

  vector<int> index_new;  // used to get the index for the splitting part
  index_new.resize(dim.size());
  for(int i=0;i<index_new.size();i++) index_new[i] = dim[i].size();
  vector< vector<int> > map;
  generate_map(map, index_new);
  index_new.resize(dim_[n].size());
  for(int i=0;i<index_new.size();i++) index_new[i] = 0;  //used to store the beg
  vector< vector<int> > sym_dim;  // for splitting part
  sym_dim.resize(map.size());
  for(int i=0;i<map.size();i++)
  {
    sym_dim[i].resize(3);
    sym_dim[i][0] = 0; // the symmetry sum
    for(int s=0;s<dim.size();s++) sym_dim[i][0] = ss(sym_dim[i][0], map[i][s]);
    sym_dim[i][1] = 1; // dimension
    for(int s=0;s<dim.size();s++) sym_dim[i][1] *= dim[s][map[i][s]];
    sym_dim[i][2] = index_new[sym_dim[i][0]]; // beg
    index_new[sym_dim[i][0]] += sym_dim[i][1];
  }
  index_new.resize(ret.dim_.size());
  tensor<T> temp;
  for(int t=0;t<sym_.size();t++) for(int s=0;s<sym_dim.size();s++)
  if(sym_dim[s][0]==sym_[t][n])
  {
    bool check = false;
    int l = 1;
    for(int i=0;i<n;i++) l *= val_[t].index_[i];
    int m = val_[t].index_[n];
    int r = 1;
    for(int i=n+1;i<dim_.size();i++) r *= val_[t].index_[i];
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
      for(int i=0;i<n;i++) index_new[i] = sym_[t][i];
      for(int i=0;i<dim.size();i++) index_new[n+i] = map[s][i];
      for(int i=n;i<dim_.size();i++) index_new[i+dim.size()] = sym_[t][i];
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
  ret.dir_ = dir_;
  ret.dim_.resize(dim_.size());
  for(int i=0;i<dim_.size();i++) ret.dim_[i] = vector<int> (1,dimension(i));
  ret.sym_.resize(1);
  ret.sym_[0] = vector<int> (dim_.size(),0);
  ret.val_.resize(1);
  ret.val_[0].index_.resize(dim_.size());
  for(int i=0;i<dim_.size();i++) ret.val_[0].index_[i] = ret.dim_[i][0];
  ret.val_[0].val_ = vector<T> (dimension(),0);
  for(int i=0;i<val_.size();i++)
  {
    vector< vector<int> > map;
    for(int j=0;j<val_[i].dimension();j++)
    {
      generate_map(map, val_[i].index_);
      int loc = 0;
      for(int s=0;s<dim_.size();s++)
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
qtensor<T> qtensor<T>::cut(int n, vector<int> cutoff) const
{
  if( n>=dim_.size() or n<0)
  {
    cerr << "Error in tensor_quantum cut: The input index number "
         << n << " is out of range 0~" << dim_.size()-1 << endl;
    return *this;
  }

  if(cutoff.size()!=dim_[n].size())
  {
    cerr << "Error in tensor_quantum cut: The input cut size "
         << cutoff.size() << " doesn't match sector number "
         << dim_[n].size() << endl;
    return *this;
  }

  qtensor<T> ret;
  ret.dir_ = dir_;
  ret.dim_.resize(dim_.size());
  for(int i=0;i<dim_.size();i++) ret.dim_[i].resize(dim_[i].size());
  ret.sym_ = sym_;
  ret.val_ = val_;

  for(int i=0;i<val_.size();i++)
  {
    ret.val_[i].resize(n, cutoff[sym_[i][n]]);
    ret.check_dim_(i);
  }
  return ret;
}


template <typename T>
qtensor<T>& qtensor<T>::plus(const qtensor& A, const qtensor& B,
                              T alpha, T beta)
{
  if( A.dim_ != B.dim_ )
  {
    cerr << "Error in tensor_quantum plus: "
            "dimensions do not match.\n" ;
    return *this;
  }
  if( A.dir_ != B.dir_ )
  {
    cerr << "Error in tensor_quantum plus: "
            "fermion directions do not match.\n" ;
    return *this;
  }
  dir_ = A.dir_;

  dim_ = A.dim_;

  vector<int> a(A.sym_.size(), 0);
  for(int i=0;i<a.size();i++) for(int j=0;j<A.dim_.size();j++)
    a[i] = a[i]*A.dim_[j].size() + A.sym_[i][j];
  vector<int> b(B.sym_.size(), 0);
  for(int i=0;i<b.size();i++) for(int j=0;j<B.dim_.size();j++)
    b[i] = b[i]*B.dim_[j].size() + B.sym_[i][j];
  vector< vector<int> > map;
  int d = 1;  // the number of all possible symmetries
  for(int i=0;i<A.dim_.size();i++) d *= A.dim_[i].size();
  for(int i=0;i<d;i++)
  {
    vector<int> loc(2,-1);
    for(int j=0;j<a.size();j++) if(a[j]==i) loc[0] = j;
    for(int j=0;j<b.size();j++) if(b[j]==i) loc[1] = j;
    if(loc[0]!=-1 or loc[1]!=-1) map.push_back(loc);
  }

  sym_.resize(map.size());
  val_.resize(map.size());
  for(int i=0;i<map.size();i++)
  {
    if(map[i][0]==-1)
    {
      sym_[i] = B.sym_[map[i][1]];
      val_[i] = B.val_[map[i][1]];
      val_[i].times(beta);
    }
    else if(map[i][1]==-1)
    {
      sym_[i] = A.sym_[map[i][0]];
      val_[i] = A.val_[map[i][0]];
      val_[i].times(alpha);
    }
    else
    {
      sym_[i] = A.sym_[map[i][0]];
      A.val_[map[i][0]].plus(val_[i], B.val_[map[i][1]], alpha, beta);
    }
  }
  return *this;
}


template <typename T>
T qtensor<T>::trace(qtensor& A, bool fermion)
{
  if( A.dim_ != dim_ )
  {
    cerr << "Error in tensor_quantum trace: "
            "dimensions do not match:\n";
    return 0.0;
  }
  if( A.dir_ != dir_ )
  {
    cerr << "Error in tensor_quantum trace: "
            "fermion directions do not match.\n" ;
    return 0.0;
  }

  T ret = 0.0;
  double sign = 1;
  for(int l=0;l<sym_.size();l++) for(int r=0;r<A.sym_.size();r++)
    if(sym_[l]==A.sym_[r])
    {
      if(fermion)
      {
        int d_b = 0;
        for(int i=0;i<A.dir_.size();i++) d_b = ss(d_b, A.sym_[r][i]*A.dir_[i]);
        int a_p = 0;
        for(int i=0;i<dir_.size();i++) if(dir_[i]==1) a_p = ss(a_p, sym_[l][i]);
        sign = 1-2*( ( ff(d_b) * ff(a_p) )%2 );
      }
      ret += val_[l].trace(A.val_[r])*sign;
    }
  return ret;
}


template <typename T>
qtensor<T>& qtensor<T>::contract(qtensor<T>& A, qtensor<T>& B,
                                 char transa, char transb, int num)
{
  int indexa = 0;
  if( 'N'==transa or 'n'==transa ) indexa = A.dim_.size()-num;
  int indexb = B.dim_.size()-num;
  if( 'N'==transb or 'n'==transb ) indexb = 0;

  // check dimension
  for(int i=0;i<num;i++) if(A.dim_[indexa+i]!=B.dim_[indexb+i])
  {
    cerr << "Error in tensor_quantum contract: "
            "The dimensions are not the same for:\n";
    cerr << "index=" << indexa+i << " A: ";
    for(int j=0;j<A.dim_[indexa+i].size();j++) cerr << A.dim_[indexa+i][j] << " ";
    cerr << endl << "index=" << indexb+i << " B: ";
    for(int j=0;j<B.dim_[indexb+i].size();j++) cerr << B.dim_[indexb+i][j] << " ";
    cerr << endl;
    return *this;
  }

  int ad = 0;
  dim_.resize(A.dim_.size()+B.dim_.size()-2*num);
  for(int i=0;i<A.dim_.size();i++) if(i<indexa or i>=indexa+num)
    dim_[ad++] = A.dim_[i];
  for(int i=0;i<B.dim_.size();i++) if(i<indexb or i>=indexb+num)
    dim_[ad++] = B.dim_[i];

  ad = 0;
  dir_.resize(A.dir_.size()+B.dir_.size()-2*num);
  for(int i=0;i<A.dir_.size();i++) if(i<indexa or i>=indexa+num)
    dir_[ad++] = A.dir_[i];
  for(int i=0;i<B.dir_.size();i++) if(i<indexb or i>=indexb+num)
    dir_[ad++] = B.dir_[i];

  sym_.clear();
  val_.clear();
  vector<int> sym(dim_.size(),0);
  vector<int> map;  // used to check if the sym already exists
  for(int i=0;i<A.val_.size();i++) for(int j=0;j<B.val_.size();j++)
  {
    bool check = true;
    for(int k=0;k<num;k++) if(A.sym_[i][indexa+k]!= B.sym_[j][indexb+k])
      check = false;

    if(check)
    {
      ad = 0;
      for(int k=0;k<A.dim_.size();k++) if(k<indexa or k>=indexa+num)
        sym[ad++] = A.sym_[i][k];
      for(int k=0;k<B.dim_.size();k++) if(k<indexb or k>=indexb+num)
        sym[ad++] = B.sym_[j][k];

      int loc = 0;
      for(int s=0;s<dim_.size();s++)
        loc = loc*dim_[s].size()+sym[s];
      bool same = false;
      for(int t=0;t<map.size();t++) if(loc==map[t])
      {
        val_[t].contract(A.val_[i], B.val_[j], transa, transb, num);
        same = true;
      }
      if(not same)
      {
        map.push_back(loc);
        sym_.push_back(sym);
        val_.resize(sym_.size());
        val_[val_.size()-1].contract(A.val_[i], B.val_[j],
                                     transa, transb, num, true);
      }
    }
  }
  for(int i=0;i<val_.size();i++) check_dim_(i);
  return *this;
}


template <typename T>
qtensor<T>& qtensor<T>::contract(const qtensor<T>& A, int a,
                                 const qtensor<T>& B, int b, bool change, bool fermion)
{
  if(A.dim_[a]!=B.dim_[b])
  {
    cerr << "Error in tensor_quantum contract: dimension doesn't match.\n"
         << "A: ";
    for(int i=0;i<A.dim_[a].size();i++) cout << A.dim_[a][i] << "\t";
    cerr << endl << "B: ";
    for(int i=0;i<B.dim_[b].size();i++) cout << B.dim_[b][i] << "\t";
    cerr << endl;
    return *this;
  }

  dim_.resize(A.dim_.size()+B.dim_.size()-2);
  for(int i=0;i<a;i++) dim_[i] = A.dim_[i];
  for(int i=a;i<a+b;i++) dim_[i] = B.dim_[i-a];
  for(int i=a+b;i<a+B.dim_.size()-1;i++) dim_[i] = B.dim_[i-a+1];
  for(int i=a+B.dim_.size()-1;i<dim_.size();i++)
    dim_[i] = A.dim_[i-B.dim_.size()+2];

  dir_.resize(A.dir_.size()+B.dir_.size()-2);
  for(int i=0;i<a;i++) dir_[i] = A.dir_[i];
  for(int i=a;i<a+b;i++) dir_[i] = B.dir_[i-a];
  for(int i=a+b;i<a+B.dir_.size()-1;i++) dir_[i] = B.dir_[i-a+1];
  for(int i=a+B.dir_.size()-1;i<dir_.size();i++)
    dir_[i] = A.dir_[i-B.dir_.size()+2];

  sym_.clear();
  val_.clear();
  vector<int> sym(dim_.size(),0);
  vector<int> map;  // used to check if the sym already exists

  for(int l=0;l<A.val_.size();l++) for(int r=0;r<B.val_.size();r++)
  if(A.sym_[l][a]==B.sym_[r][b])
  {
    for(int i=0;i<a;i++) sym[i] = A.sym_[l][i];
    for(int i=a;i<a+b;i++) sym[i] = B.sym_[r][i-a];
    for(int i=a+b;i<a+B.dim_.size()-1;i++) sym[i] = B.sym_[r][i-a+1];
    for(int i=a+B.dim_.size()-1;i<dim_.size();i++)
      sym[i] = A.sym_[l][i-B.dim_.size()+2];

    int sign = 1;
    if(fermion) if(change)
    {
      int d_a = 0;  // the sign for O_A
      for(int i=0;i<A.dir_.size();i++) d_a = ss(d_a, A.sym_[l][i]*A.dir_[i]);
      int b_p = 0;  // the sign for the right index of B
      for(int i=0;i<B.dir_.size();i++) if(B.dir_[i]==1)
        b_p = ss(b_p, B.sym_[r][i]);
      sign = (1-2*(( ff(d_a) * ff(b_p) )%2));
    }
    else
    {
      int d_b = 0;  // the sign for O_B
      for(int i=0;i<B.dir_.size();i++) d_b = ss(d_b, B.sym_[r][i]*B.dir_[i]);
      int a_p = 0;  // the sign for the right index of A
      for(int i=0;i<A.dir_.size();i++) if(A.dir_[i]==1)
        a_p = ss(a_p, A.sym_[l][i]);
      sign = (1-2*(( ff(d_b) * ff(a_p) )%2));
    }

    int loc = 0;
    for(int s=0;s<sym.size();s++) loc = loc*dim_[s].size()+sym[s];
    bool same = false;
    for(int t=0;t<map.size();t++) if(loc==map[t])
    {
      val_[t].contract(A.val_[l], a, B.val_[r], b, false, sign);
      same = true;
    }
    if(not same)
    {
      map.push_back(loc);
      sym_.push_back(sym);
      val_.resize(sym_.size());
      val_[val_.size()-1].contract(A.val_[l], a, B.val_[r], b, true, sign);
    }
  }
  for(int i=0;i<val_.size();i++) check_dim_(i);
  return *this;
}


template <typename T>
void qtensor<T>::contract(T* out, T* in,
                          const vector< vector<int> >& map,
                          char transa, char transb, int num)
{
  for(int i=0;i<map.size();i++)
    val_[map[i][0]].contract(out+map[i][2], in+map[i][1], map[i][3], map[i][4],
                     transa, transb, num, map[i][5], map[i][6]);
}


template <typename T>
bool qtensor<T>::get_map(vector< vector<int> >& map_ret,
                         vector< vector<int> >& dim_ret,
                         vector< vector<int> >& sym_ret,
                         const vector< vector<int> >& dim,
                         const vector< vector<int> >& sym,
                         char transa, char transb,
                         int num, double beta, bool back) const
{
  int indexa = 0;
  if( 'N'==transa or 'n'==transa ) indexa = dim_.size()-num;
  int indexb = dim.size()-num;
  if( 'N'==transb or 'n'==transb ) indexb = 0;

  // check dimension
  for(int i=0;i<num;i++) if(dim_[indexa+i]!=dim[indexb+i])
  {
    cerr << "Error in tensor_quantum get_map: "
            "The dimensions are not the same for i=" << i << ".\n";
    cerr << "index=" << indexa+i << " tensor: ";
    for(int j=0;j<dim_[indexa+i].size();j++) cerr << dim_[indexa+i][j] << " ";
    cerr << endl << "index=" << indexb+i << " array:  ";
    for(int j=0;j<dim[indexb+i].size();j++) cerr << dim[indexb+i][j] << " ";
    cerr << endl;
    return false;
  }

  int ad = 0;
  dim_ret.resize(dim_.size()+dim.size()-2*num);  // get dim_ret
  for(int i=0;i<dim_.size();i++) if(i<indexa or i>=indexa+num)
    dim_ret[ad++] = dim_[i];
  for(int i=0;i<dim.size();i++) if(i<indexb or i>=indexb+num)
    dim_ret[ad++] = dim[i];

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

  map_ret.clear();
  if(not back) sym_ret.clear();
  vector<int> sym_store(dim_ret.size(), 0);
  vector<int> map_store(7, 0);  // index, in_beg, out_beg, row, col, beta, alpha
  vector<int> map;  // used to check if the symmetry sector exists.
  vector<int> map_dim(1, 0);  // used to store the beginning of the sym_ret
  if(back)
  {
    map = vector<int> (sym_ret.size(), 0);
    for(int i=0;i<map.size();i++) for(int j=0;j<dim_ret.size();j++)
      map[i] = map[i]*dim_ret[j].size() + sym_ret[i][j];
    map_dim = vector<int> (sym_ret.size(), 0);
    for(int i=0;i<map_dim.size()-1;i++)
    {
      int d = 1;
      for(int j=0;j<dim_ret.size();j++) d *= dim_ret[j][sym_ret[i][j]];
      map_dim[i+1] = map_dim[i] + d;
    }
  }
  vector<bool> first(map.size(), true);  // used when back=true
  for(int i=0;i<val_.size();i++) for(int j=0;j<sym.size();j++)
  {
    bool check = true;
    for(int k=0;k<num;k++) if(sym_[i][indexa+k]!= sym[j][indexb+k])
      check = false;
    if(check and sym_dim[j][0]!=0)
    {
      map_store[0] = i;
      map_store[1] = sym_dim[j][2];
      map_store[2] = map_dim.back();  // may change
      map_store[3] = sym_dim[j][1];
      map_store[4] = sym_dim[j][0]/sym_dim[j][1];
      map_store[5] = beta;  // may change
      map_store[6] = 1;  // may change (only for fermions)

      if(not back)
      {
        int d_b = 0;
        for(int k=0;k<dir_.size();k++) d_b = ss(d_b, sym_[i][k]*dir_[k]);
        int a_p = 0;
        for(int k=0;k<dim.size()-num;k++) a_p = ss(a_p, sym[j][k]);
        map_store[6] = 1-2*(( ff(d_b) * ff(a_p) )%2);
      }

      ad = 0;
      for(int k=0;k<dim_.size();k++) if(k<indexa or k>=indexa+num)
        sym_store[ad++] = sym_[i][k];
      for(int k=0;k<dim.size();k++) if(k<indexb or k>=indexb+num)
        sym_store[ad++] = sym[j][k];
      int loc = 0;
      for(int s=0;s<dim_ret.size();s++) loc = loc*dim_ret[s].size()+sym_store[s];
      bool same = false;
      for(int t=0;t<map.size();t++) if(loc==map[t])
      {
        map_store[2] = map_dim[t];
        if(back and first[t])
          first[t] = false;
        else
          map_store[5] = 1;
        map_ret.push_back(map_store);
        same = true;
      }
      if( (not same) and (not back))
      {
        int d = 1;
        for(int s=0;s<dim_ret.size();s++)
          d *= dim_ret[s][sym_store[s]];
        d += map_dim.back();
        map_dim.push_back(d);
        map.push_back(loc);
        sym_ret.push_back(sym_store);
        map_ret.push_back(map_store);
      }
    }
  }
  return true;
}


template <typename T>
vector<double> qtensor<T>::svd(qtensor& U, qtensor<double>& S, qtensor& V,
                               int num, int cutoff)
{
  if( num>=dim_.size() or num<=0 )
  {
    cerr << "Error in tensor_quantum svd: "
            "The number of left indexes should be 1~"<< dim_.size()-1 << endl;
    return vector<double>();
  }

  qtensor<T> comb;  // combine the index to form a matrix
  comb = combine(num, dim_.size()-1);
  comb = comb.combine(0, num-1);

  int sym = ss(comb.sym_[0][0], comb.sym_[0][1]);
  bool check = true;
  for(int i=0;i<comb.sym_.size();i++)
    if(ss(comb.sym_[i][0], comb.sym_[i][1])!=sym) check = false;
  if(not check)
  {
    cerr << "Warning in tensor_quantum svd: "
            "symmetry sectors do not have constraint.\n";
    comb = comb.remove_symmetry();
    return comb.svd(U, S, V, num, cutoff);
  }

  U.dir_ = vector<int> (2, 1);
  V.dir_ = vector<int> (2, 1);
  S.dir_ = vector<int> (2, 1);

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
  if(cutoff!=0 and cutoff<sum) sum = cutoff;
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
  for(int t=sum;t<s.size();t++)  // sort the singular values beyond cutoff
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
      S.val_[i].resize(1,comb.val_[i].index_[1]);
    else S.val_[i].resize(0,comb.val_[i].index_[0]);
  }

  // if cutoff is larger than the number of non-zero singular values
  // adding the zeros for U
  if(cutoff!=0)
  {
    int diff = comb.dimension(0);  // adding the zeros for U
    diff = diff < cutoff ? diff : cutoff;  // the actual cutoff
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
      U.val_[i].resize(1, u_num[i]);
      U.dim_[1][U.sym_[i][1]] = u_num[i];
      V.val_[i].resize(0, v_num[i]);
      V.dim_[0][V.sym_[i][0]] = v_num[i];
      S.val_[i].resize(0, u_num[i]);
      S.dim_[0][S.sym_[i][0]] = u_num[i];
      S.val_[i].resize(1, v_num[i]);
      S.dim_[1][S.sym_[i][1]] = v_num[i];
    }
  }

  vector< vector<int> > dim;
  dim.resize(num);
  for(int i=0;i<num;i++) dim[i] = dim_[i];
  U = U.split(0, dim);
  dim.resize(dim_.size()-num);
  for(int i=0;i<dim.size();i++) dim[i] = dim_[i+num];
  V = V.split(1, dim);
  return s;
}


template <typename T>
int qtensor<T>::eig(double* val, T* vec, int sector)
{
  if(dim_.size()!=2)
  {
    cerr << "Error in tensor_quantum eig: cannot diagonalize. "
            "This is a " << dim_.size() << "-tensor.";
    return 0;
  }
  if(dim_[0]!=dim_[1])
  {
    cerr << "Error in tensor_quantum eig: the dimensions are not the same.\n"
            "index 0: ";
    for(int i=0;i<dim_[0].size();i++) cerr << dim_[0][i] << " ";
    cerr << endl << "index 1: ";
    for(int i=0;i<dim_[1].size();i++) cerr << dim_[1][i] << " ";
    cerr << endl;
    return 0;
  }

  for(int i=0;i<val_.size();i++) if(sym_[i][0]!=sym_[i][1])
  {
    cerr << "Warning in tensor_quantum eig: not block diagonalized."
            "Combining the symmetry sectors\n";
    qtensor<T> temp;
    temp = remove_symmetry();
    return temp.eig(val, vec, -1);
  }
  if(sector==-1)
  {
    if(sym_.size()==1) return val_[0].eig(val, vec);
    qtensor<T> temp;
    temp = remove_symmetry();
    return temp.eig(val, vec, -1);
  }
  int loc = -1;
  for(int i=0;i<sym_.size();i++) if(sym_[i][0]==sector) loc = i;
  if(loc==-1)
  {
    cerr << "Warning in tensor_quantum eig: cannot find symmetry sector "
         << sector << endl;
    int d = dim_[0][sector];
    for(int i=0;i<d;i++) val[i] = 0;
    for(int i=0;i<d;i++) vec[i] = 0;
    vec[0] = 1;
    return 0;
  }
  return val_[loc].eig(val, vec);
}


template class qtensor<double>;
template class qtensor< complex<double> >;

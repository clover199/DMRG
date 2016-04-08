
#include "global.h"
#include "mps.h"
#include "mpo.h"

vector<string> set_para_name(string s1, string s2, string s3, string s4,
                             string s5, string s6, string s7, string s8)
{
  vector<string> ret;
  ret.clear();
  if(s1!="") { ret.push_back(s1);
  if(s2!="") { ret.push_back(s2);
  if(s3!="") { ret.push_back(s3);
  if(s4!="") { ret.push_back(s4);
  if(s5!="") { ret.push_back(s5);
  if(s6!="") { ret.push_back(s6);
  if(s7!="") { ret.push_back(s7);
  if(s8!="") { ret.push_back(s8);}}}}}}}}
  return ret;
}


// This function is only for use in this file
bool split(string & str, string sep_char, string & str_bef, string & str_aft)
{
  char buffer[50000];

  size_t found = str.find(sep_char);
  if (found != string::npos)
  {
    int loc = int(found);

    size_t length = str.copy(buffer,loc,0);
    buffer[length]='\0';
    str_bef=buffer;

    length = str.copy(buffer,str.size()-loc,loc+1);
    buffer[length]='\0';
    str_aft=buffer;
    return true;
  }
  else return false;
}


vector<double> set_para_val(int argc_, char *argv_[], int& sites, int& cutoff, int& sweep,
                            int& sector, vector<string>& para_name, string& filename)
{
  vector<double> ret(para_name.size(),0);
  if(ret.size()) ret[0] = 1;
  string str_bef, str_aft;
  for(int i=0;i<argc_;i++)
  {
    string str_input(argv_[i]);
    split(str_input, "=", str_bef, str_aft);
    if(str_bef.compare("l") == 0)
      sites = (int)atof(str_aft.c_str());
    if(str_bef.compare("m") == 0)
      cutoff = (int)atof(str_aft.c_str());
    if(str_bef.compare("n") == 0)
      sweep = (int)atof(str_aft.c_str());
    if(str_bef.compare("s") == 0)
      sector = (int)atof(str_aft.c_str());
    for(int j=0;j<para_name.size();j++) if(str_bef.compare(para_name[j]) == 0)
      ret[j] = atof(str_aft.c_str());
  }
  
  cout << "\n ***** The parameters for this running ***** \n\n"
          "  sites(l) = " << sites << "\t cutoff(m) = " << cutoff
       << "\t sweep(n) = " << sweep << "\t symmetry_section(s) = " << sector << endl;
  for(int i=0;i<ret.size();i++) cout << "\t " << para_name[i] << " = " << ret[i] << endl;
  cout << endl;
  
  std::stringstream ss;
  ss << "_L_" << sites << "_m_" << cutoff << "_n_" << sweep << "_s_" << sector;
  for(int i=0;i<ret.size();i++) ss << "_" << para_name[i] << "_" << ret[i];
  filename = ss.str()+"_.txt";
  
  return ret;
}


template <typename T>
qtensor<T> combine_mpo(qtensor<T>& H1, qtensor<T>& H2)
{
  qtensor<T> ret;
  ret.contract(H1, 2, H2, 0);
  ret = ret.exchange(4,5);
  ret = ret.combine(4,5);
  ret = ret.combine(1,2);
  return ret;
}
template qtensor<double> combine_mpo(qtensor<double>& H1, qtensor<double>& H2);
template qtensor<complex<double> > combine_mpo(qtensor<complex<double> >& H1, qtensor<complex<double> >& H2);



template <typename T>
void two_sites(int l, int r, mps<T>& my_mps, int sector, 
               ofstream& data_energy, ofstream& data_singular);
template <typename T>
void grow2right(int l, mps<T>& my_mps, mpo<T>& my_mpo);
template <typename T>
void grow2left(int r, mps<T>& my_mps, mpo<T>& my_mpo);


template <typename T>
void ed(mps<T>& my_mps, mpo<T>& my_mpo, int cut, int swp, int sector, const string& filename)
{
  int L = my_mps.size();

  ofstream data_energy, data_singular;
  string name = "energy"+filename;
  data_energy.open(name.c_str());
  name = "singular"+filename;
  data_singular.open(name.c_str());

  
  for(int i=0;i<L/2;i++)
  {
    grow2right(i, my_mps, my_mpo);
    grow2left(L-i-1, my_mps, my_mpo);
  }
  if(L%2)
    grow2left(L/2, my_mps, my_mpo);
  cout << "********** Calculating (L=" << L << ")... **********\n";
  two_sites(L/2-1, L/2, my_mps, sector, data_energy, data_singular);

  data_energy.close();
  data_singular.close();
}
template void ed(mps<double>& my_mps, mpo<double>& my_mpo,
                 int cut, int swp, int sector, const string& filename);
template void ed(mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo,
                 int cut, int swp, int sector, const string& filename);



template <typename T>
double two_sites(int l, int r, int cutoff,
               mps<T>& my_mps, mpo<T>& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular);
template <typename T>
double update_two(int l, int r, int cutoff,
                mps<T>& my_mps, mpo<T>& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular);
template <typename T>
double move2right(int l, int r, int cutoff,
                mps<T>& my_mps, mpo<T>& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular);
template <typename T>
double move2left(int l, int r, int cutoff,
               mps<T>& my_mps, mpo<T>& my_mpo,
               int sector, ofstream& data_energy, ofstream& data_singular);


template <typename T>
void dmrg(mps<T>& my_mps, mpo<T>& my_mpo,
          int cutoff, int sweep, int sector, const string& filename)
{
  int L = my_mps.size();

  ofstream data_energy, data_singular;
  string name = "energy"+filename;
  data_energy.open(name.c_str());
  name = "singular"+filename;
  data_singular.open(name.c_str());

  int pre_cutoff = 10;
  if(sweep==-1) pre_cutoff = cutoff;
  int pre_sweep = cutoff/100;
  
  cout << "********** starting l=" << 0 << " r=" << L-1 << " **********\n";
  two_sites(0, L-1, cutoff*10, my_mps, my_mpo, sector, data_energy, data_singular);
  for(int i=1;i<L/2;i++)
  {
    cout << ">>>>>>>>>>> starting l=" << i << " r=" << L-1-i << " <<<<<<<<<<\n";
    update_two(i, L-i-1, pre_cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
  }
  if(sweep!=-1) for(int i=L/2;i<L-1;i++)
  {
    cout << ">>>>>>>>>> starting l=" << i << " r=" << i+1 << " >>>>>>>>>>\n";
    move2right(i, i+1, pre_cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
  }
  if(sweep!=-1) for(int s=0;s<pre_sweep-1;s++)
  {
    for(int i=L-2;i>0;i--)
    {
      cout << "<<<<<<<<<< pre-sweep: cutoff=" << 100*(s+1)
           << " l=" << i-1 << " r=" << i << " <<<<<<<<<<\n";
      move2left(i-1, i, 100*(s+1), my_mps, my_mpo, sector, data_energy, data_singular);
    }
    for(int i=1;i<L-1;i++)
    {
      cout << ">>>>>>>>>> pre-sweep: cutoff=" << 100*(s+1)
           << " l=" << i << " r=" << i+1 << " >>>>>>>>>>\n";
      move2right(i, i+1, 100*(s+1), my_mps, my_mpo, sector, data_energy, data_singular);
    }
  }
  for(int s=0;s<sweep;s++)
  {
    for(int i=L-2;i>0;i--)
    {
      cout << "<<<<<<<<<< sweep=" << s+1
           << " l=" << i-1 << " r=" << i << " <<<<<<<<<<\n";
      move2left(i-1, i, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
    }
    for(int i=1;i<L-1;i++)
    {
      cout << ">>>>>>>>>> sweep=" << s+1
           << " l=" << i << " r=" << i+1 << " >>>>>>>>>>\n";
      move2right(i, i+1, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
    }
  }
  if(sweep!=-1) for(int i=L-2;i>=L/2;i--)
  {
    cout << "<<<<<<<<<< final l=" << i-1 << " r=" << i << " <<<<<<<<<<\n";
    move2left(i-1, i, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
  }
  // I don't know why, but only adding this can make the MPS orthogonal
//  move2right(L/2-1, L/2, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);


  data_energy.close();
  data_singular.close();
}
template void dmrg(mps<double>& my_mps, mpo<double>& my_mpo,
          int cutoff, int sweep, int sector, const string& filename);
template void dmrg(mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo,
          int cutoff, int sweep, int sector, const string& filename);


template <typename T>
void dmrg2(mps<T>& my_mps, mpo<T>& my_mpo,
           int cutoff, int sweep, int sector, const string& filename)
{
  int L = my_mps.size();

  ofstream data_energy, data_singular;
  string name = "energy"+filename;
  data_energy.open(name.c_str());
  name = "singular"+filename;
  data_singular.open(name.c_str());

  int pre_cutoff = 10;
  if(sweep==-1) pre_cutoff = cutoff;
  int pre_sweep = cutoff/100;

  cout << "********** starting l=" << 0 << " r=" << L-1 << " **********\n";
  two_sites(0, L-1, 0, my_mps, my_mpo, sector, data_energy, data_singular);
  for(int i=1;i<L/2;i++)
  {
    cout << ">>>>>>>>>>> starting l=" << i << " r=" << L-1-i << " <<<<<<<<<<\n";
    update_two(i, L-i-1, pre_cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
  }
  if(sweep!=-1) for(int i=L/2;i<L-2;i++)
  {
    cout << ">>>>>>>>>> starting l=" << i << " r=" << i+1 << " >>>>>>>>>>\n";
    update_two(i, i+1, pre_cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
  }
  if(sweep!=-1) for(int s=0;s<pre_sweep-1;s++)
  {
    for(int i=L-3;i>1;i--)
    {
      cout << "<<<<<<<<<< pre-sweep: cutoff=" << 100*(s+1)
           << " l=" << i-1 << " r=" << i << " <<<<<<<<<<\n";
      update_two(i-1, i, 100*(s+1), my_mps, my_mpo, sector, data_energy, data_singular);
    }
    for(int i=2;i<L-2;i++)
    {
      cout << ">>>>>>>>>> pre-sweep: cutoff=" << 100*(s+1)
           << " l=" << i << " r=" << i+1 << " >>>>>>>>>>\n";
      update_two(i, i+1, 100*(s+1), my_mps, my_mpo, sector, data_energy, data_singular);
    }
  }
  for(int s=0;s<sweep;s++)
  {
    for(int i=L-3;i>1;i--)
    {
      cout << "<<<<<<<<<< sweep=" << s+1
           << " l=" << i-1 << " r=" << i << " <<<<<<<<<<\n";
      update_two(i-1, i, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
    }
    for(int i=2;i<L-2;i++)
    {
      cout << ">>>>>>>>>> sweep=" << s+1
           << " l=" << i << " r=" << i+1 << " >>>>>>>>>>\n";
      update_two(i, i+1, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
    }
  }
  if(sweep!=-1) for(int i=L-3;i>=L/2;i--)
  {
    cout << "<<<<<<<<<< final l=" << i-1 << " r=" << i << " <<<<<<<<<<\n";
    update_two(i-1, i, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
  }
  data_energy.close();
  data_singular.close();
}
template void dmrg2(mps<double>& my_mps, mpo<double>& my_mpo,
          int cutoff, int sweep, int sector, const string& filename);
template void dmrg2(mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo,
          int cutoff, int sweep, int sector, const string& filename);




template <typename T>
void idmrg(mps<T>& my_mps, mpo<T>& my_mpo,
          int cutoff, int sweep, int sector, const string& filename)
{
  int L = my_mps.size();
  double ene = 0;

  ofstream data_energy, data_singular;
  string name = "energy"+filename;
  data_energy.open(name.c_str());
  name = "singular"+filename;
  data_singular.open(name.c_str());
  
  cout << "********** L=" << 2 << " **********\n";
  ene = two_sites(0, L-1, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
  ene = ene/2;
  int l = 1;
  double diff = 1;
  while(l<L/2 and abs(diff)>TOL*100)
  {
    double store = 0;
    cout << ">>>>>>>>>>> L=" << l*2+2 << " <<<<<<<<<<\n";
    store = update_two(l, L-l-1, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
    store = store / (l*2+2);
    my_mps.clear(l-1);
    my_mps.clear(L-l);
    diff = ene - store;
    cout << "Energy difference: " << diff << endl;
    ene = store;
    l++;
  }

  data_energy.close();
  data_singular.close();
}
template void idmrg(mps<double>& my_mps, mpo<double>& my_mpo,
          int cutoff, int sweep, int sector, const string& filename);
template void idmrg(mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo,
          int cutoff, int sweep, int sector, const string& filename);




template <typename T>
void idmrg1(mps<T>& my_mps, mpo<T>& my_mpo,
          int cutoff, int sweep, int sector, const string& filename)
{
  int L = my_mps.size();
  double ene = 0;

  ofstream data_energy, data_singular;
  string name = "energy"+filename;
  data_energy.open(name.c_str());
  name = "singular"+filename;
  data_singular.open(name.c_str());
  
  cout << "********** L=" << 2 << " **********\n";
  ene = two_sites(0, L-1, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
  int l = 1;
  double diff = 1;
  while(l<L/2 and abs(diff)>TOL)
  {
    double store = 0;
    cout << ">>>>>>>>>>> L=" << l*2+1 << " <<<<<<<<<<\n";
    move2left(l-1, L-l-1, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
    my_mps.clear(L-l);
    cout << ">>>>>>>>>>> L=" << l*2+2 << " <<<<<<<<<<\n";
    store = move2right(l, L-l-1, cutoff, my_mps, my_mpo, sector, data_energy, data_singular);
    store = store / (l*2+2);
    my_mps.clear(l-1);
    diff = ene - store;
    cout << "Energy difference: " << diff << endl;
    ene = store;
    l++;
  }

  data_energy.close();
  data_singular.close();
}
template void idmrg1(mps<double>& my_mps, mpo<double>& my_mpo,
          int cutoff, int sweep, int sector, const string& filename);
template void idmrg1(mps<complex<double> >& my_mps, mpo<complex<double> >& my_mpo,
          int cutoff, int sweep, int sector, const string& filename);



template <typename T>
T calc(mps<T>& my_mps)
{
  int sites = my_mps.size();
  if(sites<2)
  {
    cerr << "Error in calc(mps): empty MPS.\n";
    return 0;
  }

  qtensor<T> left;
  qtensor<T> end;
  qtensor<T> temp;

  temp = my_mps[0];
  temp.conjugate();
  left.contract(temp, my_mps[0], 'T', 'N', 1);
  for(int i=1;i<sites-1;i++)
  {
    temp = my_mps[i];
    end.contract(left, temp, 'N', 'N', 1);
    temp.conjugate();
    left.contract(temp, end, 'T', 'N', 2);
  }
  temp = my_mps[sites-1];
  temp.conjugate();
  end.contract(temp, my_mps[sites-1], 'N', 'T', 1);
  return left.trace(end);
}

template double calc(mps<double>& my_mps);
template complex<double> calc(mps< complex<double> >& my_mps);


template <typename T>
T calc(mps<T>& l_mps, mps<T>& r_mps)
{
  if(l_mps.size()!=r_mps.size())
  {
    cout << "Error in calc(mps, mps): sizes do not match: "
         << l_mps.size() << "!=" << r_mps.size() << endl;
    return 0;
  }
  int sites = l_mps.size();
  if(sites<2)
  {
    cerr << "Error in calc(mps, mps): empty MPS.\n";
    return 0;
  }

  qtensor<T> left;
  qtensor<T> temp;
  qtensor<T> end;

  temp = r_mps[0];
  temp.conjugate();
  left.contract(temp, l_mps[0], 'T', 'N', 1);
  for(int i=1;i<sites-1;i++)
  {
    end.contract(left, l_mps[i], 'N', 'N', 1);
    temp = r_mps[i];
    temp.conjugate();
    left.contract(temp, end, 'T', 'N', 2);
  }
  temp = r_mps[sites-1];
  temp.conjugate();
  end.contract(temp, l_mps[sites-1], 'N', 'T', 1);
  return left.trace(end);
}

template double calc(mps<double>& l_mps, mps<double>& r_mps);
template complex<double> calc(mps< complex<double> >& l_mps,
                              mps< complex<double> >& r_mps);


template <typename T>
T calc(mps<T>& my_mps, mpo<T>& my_mpo)
{
  int sites = my_mps.size();
  if(sites!=my_mpo.size())
  {
    cout << "Error in calc(mps, mpo): sizes do not match: "
         << sites << "!=" << my_mpo.size() << endl;
    return 0;
  }
  if(sites<2)
  {
    cerr << "Error in calc(mps, mpo): empty MPS.\n";
    return 0;
  }

  qtensor<T> left;
  qtensor<T> end;
  qtensor<T> temp;

  end.contract(my_mpo[0], my_mps[0], 'N', 'N', 1);
  temp = my_mps[0];
  temp.conjugate();
  left.contract(temp, end, 'T', 'N', 1);
  for(int i=1;i<sites-1;i++)
  {
    end.contract(left, 1, my_mpo[i], 0);
    left = end.exchange(3,4);
    end.contract(left, my_mps[i], 'N', 'N', 2);
    temp = my_mps[i];
    temp.conjugate();
    left.contract(temp, end, 'T', 'N', 2);
  }
  qtensor<T> store;
  store.contract(my_mpo[sites-1], my_mps[sites-1], 'N', 'T', 1);
  temp = my_mps[sites-1];
  temp.conjugate();
  end.contract(temp, store, 'N', 'N', 1);
  return left.trace(end, true);
}

template double calc(mps<double>& my_mps, mpo<double>& my_mpo);
template complex<double> calc(mps< complex<double> >& my_mps,
                              mpo< complex<double> >& my_mpo);

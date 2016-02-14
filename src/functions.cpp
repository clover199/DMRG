
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


// only used inside this file
qtensor<double> connect(qtensor<double>& cent, qtensor<double>& right)
{
  qtensor<double> ret;
  ret = cent * right;
  return ret;
}


// only used inside this file
qtensor<complex<double> > connect(qtensor<double>& cent, qtensor<complex<double> >& right)
{
  qtensor<complex<double> > ret;
  ret = cent.comp() * right;
  return ret;
}


// only used inside this file
qtensor<complex<double> > connect(qtensor<complex<double> >& left, qtensor<double>& cent)
{
  qtensor<complex<double> > ret = cent.comp();
  ret = left * ret;
  return ret;
}


template <typename T>
T calc(mps<T>& my_mps)
{
  int loc = my_mps.position();
  int sites = my_mps.size();
  if(sites<0 or loc<0) return 0;

  qtensor<T> left;
  qtensor<T> right;
  qtensor<T> temp;
  qtensor<T> store;

  temp = my_mps[0];
  temp.conjugate();
  left.contract(temp, my_mps[0], 'T', 'N', 1);
  for(int i=1;i<loc;i++)
  {
    store.contract(left, my_mps[i], 'N', 'N', 1);
    temp = my_mps[i];
    temp.conjugate();
    left.contract(temp, store, 'T', 'N', 2);
  }
  temp = my_mps[sites-1];
  temp.conjugate();
  right.contract(temp, my_mps[sites-1], 'N', 'T', 1);
  for(int i=sites-2;i>loc;i--)
  {
    temp = my_mps[i];
    temp.conjugate();
    store.contract(temp, right, 'N', 'N', 1);
    right.contract(store, my_mps[i], 'N', 'T', 2);
  }
  qtensor<double> singular = my_mps.center(loc);
  temp = connect(singular, my_mps[loc]);
  temp.conjugate();
  store.contract(temp, right, 'N', 'N', 1);
  temp.conjugate();
  right.contract(store, temp, 'N', 'T', 2);
  return left.trace(right);
}

template double calc(mps<double>& my_mps);
template complex<double> calc(mps< complex<double> >& my_mps);


template <typename T>
T calc(mps<T>& l_mps, mps<T>& r_mps)
{
  if(l_mps.size()<=0 or r_mps.size()<=0) return 0;
  if(l_mps.position()<=0 or r_mps.position()<=0) return 0;
  if(l_mps.size()!=r_mps.size())
  {
    cout << "Error in calc(mps, mps): sizes do not match: "
         << l_mps.size() << "!=" << r_mps.size() << endl;
    return 0;
  }

  qtensor<T> left;
  qtensor<T> temp;
  qtensor<T> store;

  temp = r_mps[0];
  temp.conjugate();
  left.contract(temp, l_mps[0], 'T', 'N', 1);
  for(int i=1;i<l_mps.size()-1;i++)
  {
    if(i==l_mps.position()-1) temp = connect(l_mps[i], l_mps.center());
    else temp = l_mps[i];
    store.contract(left, temp, 'N', 'N', 1);
    temp = r_mps[i];
    if(i==r_mps.position()-1) temp = connect(r_mps[i], r_mps.center());
    else temp = r_mps[i];
    temp.conjugate();
    left.contract(temp, store, 'T', 'N', 2);
  }
  temp = r_mps[r_mps.size()-1];
  temp.conjugate();
  store.contract(temp, l_mps[l_mps.size()-1], 'N', 'N', 1);
  return left.trace(store);
}

template double calc(mps<double>& l_mps, mps<double>& r_mps);
template complex<double> calc(mps< complex<double> >& l_mps,
                              mps< complex<double> >& r_mps);


template <typename T>
T calc(mps<T>& my_mps, mpo<T>& my_mpo)
{
  int loc = my_mps.position();
  int sites = my_mps.size();
  if(sites<0 or loc<0) return 0;

  qtensor<T> left;
  qtensor<T> right;
  qtensor<T> temp;
  qtensor<T> store;

  store.contract(my_mpo[0], my_mps[0], 'N', 'N', 1);

  temp = my_mps[0];
  temp.conjugate();
  left.contract(temp, store, 'T', 'N', 1);
  for(int i=1;i<loc;i++)
  {
    store.contract(left, 1, my_mpo[i], 0);
    left = store.exchange(3,4);
    store.contract(left, my_mps[i], 'N', 'N', 2);
    temp = my_mps[i];
    temp.conjugate();
    left.contract(temp, store, 'T', 'N', 2);
  }
  store.contract(my_mpo[sites-1], my_mps[sites-1], 'N', 'T', 1);
  temp = my_mps[sites-1];
  temp.conjugate();
  right.contract(temp, store, 'N', 'N', 1);
  for(int i=sites-2;i>loc;i--)
  {
    temp = my_mpo[i].exchange(0,2);
    store.contract(right, 1, temp, 0, true);
    right = store.exchange(0,1);
    temp = my_mps[i];
    temp.conjugate();
    store.contract(temp, right, 'N', 'N', 2);
    right.contract(store, my_mps[i], 'N', 'T', 2);
  }
  temp = my_mpo[loc].exchange(0,2);
  store.contract(right, 1, temp, 0, true);
  right = store.exchange(0,1);
  qtensor<double> singular = my_mps.center(loc);
  temp = connect(singular, my_mps[loc]);
  temp.conjugate();
  store.contract(temp, right, 'N', 'N', 2);
  temp.conjugate();
  right.contract(store, temp, 'N', 'T', 2);

  return left.trace(right, true);
}

template double calc(mps<double>& my_mps, mpo<double>& my_mpo);
template complex<double> calc(mps< complex<double> >& my_mps,
                              mpo< complex<double> >& my_mpo);

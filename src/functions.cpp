
#include "global.h"

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
                            vector<string>& para_name)
{
  vector<double> ret(para_name.size(),1);
  string str_bef, str_aft;
  for(int i=0;i<argc_;i++)
  {
    string str_input(argv_[i]);
    split(str_input, "=", str_bef, str_aft);
    if(str_bef.compare("L") == 0)
      sites = (int)atof(str_aft.c_str());
    if(str_bef.compare("m") == 0)
      cutoff = (int)atof(str_aft.c_str());
    if(str_bef.compare("n") == 0)
      sweep = (int)atof(str_aft.c_str());
    for(int j=0;j<para_name.size();j++) if(str_bef.compare(para_name[j]) == 0)
      ret[j] = atof(str_aft.c_str());
  }
  
  cout << "\n ***** The parameters for this running ***** \n\n"
          "  sites(L) = " << sites << "\t cutoff(m) = " << cutoff << "\t sweep(n) = " << sweep << endl;
  for(int i=0;i<ret.size();i++) cout << "\t " << para_name[i] << " = " << ret[i] << endl;
  cout << endl;
  
  std::stringstream ss;
  ss << "_L_" << sites << "_m_" << cutoff << "_n_" << sweep;
  for(int i=0;i<ret.size();i++) ss << "_" << para_name[i] << "_" << ret[i];
  filename = ss.str()+"_.txt";
  
  return ret;
}

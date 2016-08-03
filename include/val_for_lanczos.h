#ifndef _VAL_FOR_LANCZOS_
#define _VAL_FOR_LANCZOS_

#include "global.h"
#include "qtensor.h"

template <typename T>
class lanczos
{
public:
  qtensor<T> lenv;  // the left environment
  qtensor<T> renv;  // the right environment
  T * store;
  vector< vector<int> > lmap;
  vector< vector<int> > rmap;
  int r_num;
  int l_num;
#ifdef PBC
  qtensor<T> ledge;  // the left edge
  qtensor<T> redge;  // the right edge
  vector< vector<int> > lemap;
  vector< vector<int> > remap;
#endif
};

#endif
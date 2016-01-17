#include "global.h"

// generate the full array of all possible symmetries (sorted) given index 
// i.e. index=[4,1,3], then
// map=[ [0,0,0]
//       [0,0,1]
//       [0,0,2]
//       [1,0,0]
//       [1,0,1]
//       [1,0,2]
//       [2,0,0]
//       [2,0,1]
//       [2,0,2]
//       [3,0,0]
//       [3,0,1]
//       [3,0,2] ]
void generate_map(vector< vector<int> >& map, const vector<int>& index)
{
  int dim = 1;  // dimension of map
  for(int i=0;i<index.size();i++) dim *= index[i];
  map.resize(dim);
  for(int i=0;i<dim;i++) map[i].resize(index.size());
  int count = 1;  // the periodicity of the symmetry for the index i
  for(int i=0;i<index.size();i++)
  {
    for(int p=0;p<count;p++)
      for(int s=0;s<index[i];s++)
        for(int k=0;k<dim/count/index[i];k++)
          map[k+s*dim/count/index[i]+p*dim/count][i] = s;
    count *= index[i];
  }
}
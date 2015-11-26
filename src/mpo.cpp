
#include "mpo.h"

mpo::mpo(int sites, tensor ham)
{
  sites_ = sites;
  A_.push_back(ham);
}

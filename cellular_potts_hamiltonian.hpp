#ifndef __HAMILTONIAN__
#define __HAMILTONIAN__
#include <vector>
#include "cellular_potts_definition.hpp"
#include "cellular_potts_type.hpp"
#include "cellular_potts_site.hpp"
class hamiltonian_system_class
{
};
class geometry_cellular_potts
{
  /*======================
    Members
   =======================*/
private: std::vector<long long int> neighbor_sites;
  /*======================
    Methods
   =======================*/
  geometry_cellular_potts();
  /*======================
    Constructor
   =======================*/
};
#endif // __HAMILTONIAN__

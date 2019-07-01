#ifndef __FIELD__
#define __FIELD__
#include <vector>
#include "toolbox.hpp"
#include "cellular_potts_definition.hpp"
// don't include cellular_potts_cell.hpp for saving dependency
#include "cellular_potts_io.hpp"
class region_cellular_potts_class
{
  /*======================
    Members
   =======================*/
private: std::string region_identifier;
private: std::vector<long int> lower_boundary;
private: std::vector<long int> upper_boundary;
private: std::vector<double> external_field;
  /*======================
    Methods
   =======================*/
private: void set_region_identifier(
				    const std:: string & region_identifier
				    );
private: std::string get_region_identifier() const;
  //
private: void set_boudary(
			  const std::vector<long int> & input_boundary
			  );
private: long int get_bounary(
			      const std::string & bound_identifier,
			      const int & direction_index
			      ) const ;
  //
private: void set_external_field(
				 const std::vector<long int> & field
				 );
private: double get_field(
			  int & direction_index
			  )const;
  /*======================
    Constructor
   =======================*/
}
#endif __FIELD__

#ifndef __REGION__
#define __REGION__
#include <vector>
#include "toolbox.hpp"
#include "cellular_potts_definition.hpp"
// don't include cellular_potts_cell.hpp for saving dependency
#include "cellular_potts_io.hpp"
#include "cellular_potts_site.hpp"
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
public: void set_region_identifier(
				   const std::string & input_region_identifier
				   );
public: std::string get_region_identifier() const;
  //
public: void set_boundary(
			  const std::string & boundary_identifier,
			  const std::vector<long int> & input_boundary
			  );
public: long int get_boundary(
			      const std::string & boundary_identifier,
			      const int & direction_index
			      ) const ;
  //
public: void set_external_field(
				const std::vector<double> & input_field
				);
public: double get_field(
			  int & direction_index
			  ) const;
  //
public: void show() const;
  /*======================
    Constructor
   =======================*/
public: region_cellular_potts_class();  
};
//
class region_system_class{
  /*======================
    Members
    =======================*/
private: std::vector<region_cellular_potts_class> regions;
private: io_cellular_potts io_method;
private: std::string structure_item;
private: long int space_dimension;
private: int number_of_regions;
  // work_memory
private: std::vector<long int> work_vector_for_coordinate;// size=space_dimension
private: int in_region_flag;
private: int out_region_flag;
private: int error_region_flag;
  /*======================
    Methods
   =======================*/
private: void initialize_region();
  //
public: void get_field(
		       const site_system_class & site_system,
		       const long long int & site_index,
		       std::vector<double> & field
		       ) ;
  //
private: void set_structure_region_item(
					const int & type_index,
					const std::string & child
					);
  //
private: int region_identify(
			     const int region_index,
			     const site_system_class & site_system,
			     const long long int & site_index
			     );
private: void show_regions() ;
  //
public: void show_region(
			 const int & region_index
			 ) ;
  //
public: void get_region_field(
			      const int & region_index,
			      std::vector<double> & field
			      ) const;
  //
public: void set_region_field(
			      const int & region_index,
			      const std::vector<double> & field
			      );
  //
public: void test_region(
			  const site_system_class & site_system,
			  const long long int & site_index
			  );
  /*======================
    Constructor
    =======================*/
public: region_system_class(
			    const model_parameters_cellular_potts_class & model
			    );
};
#endif // __REGION__

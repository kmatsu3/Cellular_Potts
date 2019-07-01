#ifndef __ADHESION_TYPE__
#define __ADHESION_TYPE__
#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
#include "cellular_potts_type.hpp"
#include "toolbox.hpp"
#include <vector>
//
class adhesion_binding_table{
  /*======================
    Members
    =======================*/
private: int adhesion_index;
private: std::vector<long int> bind_partner_index;
private: std::vector<long int> bind_partner_pointer;
public: std::vector<double> bind_coupling;
private: io_cellular_potts io_method;
private: toolbox tools;
  /*======================
    Methods
   =======================*/
public: void set_table(
		       const int & input_adhesion_index,
		       const long int & number_of_cells,
		       const double & coupling_constant
		       );
public: adhesion_binding_table();
private: void show_binding_table(
				 const long int & number_of_cells
				 );
private: std::string set_structure_bind_table_item(
						   const int & adhesion_index,
						   const long int &cell_index,
						   const std::string & child
						   );
};
//
class adhesion_cellular_potts_class{
  /*======================
    Members
   =======================*/
private: int cell_type_1;
private: int polarity_type_1;
private: std::vector<double> adhesion_polarity_1;
private: int cell_type_2;
private: int polarity_type_2;
private: std::vector<double> adhesion_polarity_2;
private: double coupling_constant;
private: std::string interaction_type;
public: adhesion_binding_table bind_table;
  /*======================
    Methods
   =======================*/
public: void set_adhesion_int_value(
				    const std::string & type,
				    const int & value
				    );
public: void set_adhesion_double_value(
				       const std::string & type,
				       const double & value
				       );
public: void set_adhesion_double_vector(
					const int & input_cell_index,
					const std::vector<double> & adhesion_polarity
					);
public: void set_adhesion_string(
				 const std::string & type,
				 const std::string & value
				 );
public: void get_adhesion_polarity(
				   const int & input_cell_index,
				   std::vector<double> & polarity
				   ) const ;
public: double get_adhesion_cell1_component(
					    const int & component_index
					    ) const ;
public: double get_adhesion_cell2_component(
					    const int & component_index
					    ) const ;
public: int get_adhesion_cell_type(
				   const int & input_cell_index
				   ) const ;
public: int get_adhesion_polarity_type(
				       const int & input_cell_index
				       ) const ;
public: double get_coupling_constant() const;
public: std::string get_interaction_type() const;
  /*======================
    Constructor
   =======================*/
public: adhesion_cellular_potts_class();
};
//
class adhesion_system_class{
  /*======================
    Members
   =======================*/
  /* model parameters */
private: int space_dimension;
public: long int number_of_cells;
public: int number_of_cell_types;
public: int number_of_adhesion;
public: int number_of_adhesion_components;
private: std::vector<adhesion_cellular_potts_class> adhesions;
public: std::vector<bool> typepair_to_adhesion_flags;
public: std::vector<int> typepair_to_adhesion_pointers;
private: std::vector<int> typepair_to_adhesion_maps;
private: std::vector<int> type_and_type_adhesion_pointer;
public: std::vector<double> coupling_constants;
public: std::vector<double> polar_components_1; 
public: std::vector<double> polar_components_2; 
private: io_cellular_potts io_method;
  // adhesion_keywords
private: std::string adhesion_type_normal;
private: std::string adhesion_type_tight;
  /*======================
    Methods
   =======================*/
public: void initialize_adhesion();
private: std::string set_structure_adhesion_item(
						 const int & adhesion_index,
						 const std::string & child
						 );
public: void get_adhesion_polarity(
				   const int & adhesion_index,
				   const int & cell_index,
				   std::vector<double> & polarity
				   ) const ;
public: double get_adhesion_cell1_component(
					    const int & adhesion_index,
					    const int & component_index
					    ) const;
public: double get_adhesion_cell2_component(
					    const int & adhesion_index,
					    const int & component_index
					    ) const;
public: double get_coupling_constant(
				     const int & adhesion_index
				     ) const ;
public: void set_coupling_constant(
				   const int & adhesion_index,
				   const double & value
				   );
public: std::string get_interaction_type(
					 const int & adhesion_index
					 ) const ;
public: int get_adhesion_cell_type(
				   const int & adhesion_index,
				   const int & cell_index
				   ) const ;
public: int get_adhesion_polarity_type(
				       const int & adhesion_index,
				       const int & cell_index
				       ) const ;
private: void interaction_key_check(const std::string & value);
public: void make_typepair_to_adhesion_maps();
public: void show_adhesion();
public: void show_map_table();
  /*======================
    Constructor
   =======================*/
public : adhesion_system_class(
			       const model_parameters_cellular_potts_class & model,
			       const type_system_class & type_system
			       );
};
//
#endif //__ADHESION_TYPE__

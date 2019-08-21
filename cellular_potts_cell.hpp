#ifndef __CELL__
#define __CELL__
#include  "cellular_potts_site.hpp"
#include  "cellular_potts_definition.hpp"
#include  "cellular_potts_io.hpp"
#include  "cellular_potts_type.hpp"
#include  "toolbox.hpp"
#include <vector>
//
class cell_cellular_potts_class {
  /*======================
    Members
   =======================*/
private: int Type;
  /*======================
    Methods
   =======================*/
  //
public: void set_type(const int & input_type);
public: int get_type() const;
  /*======================
    Constructor
   =======================*/
public:cell_cellular_potts_class();
};
class cell_system_class{
  /*======================
    Members
   =======================*/
private: std::vector<cell_cellular_potts_class> cells;
private: long int buffer_cell;
private: int int_error_number;
  /* model parameters */
private: int space_dimension;
private: long long int number_of_sites;
private: long int number_of_cells;
private: int number_of_cell_types;
private: std::vector<long int> system_dimensions;
public: std::vector<bool> mobile_table;
  /* Work Memory*/
private: std::vector<long int> work_positions;
  /* Subclass */
private: toolbox tool;
private: io_cellular_potts io_method;
  /*======================
    Methods
   =======================*/
private: int type_define(
			 const model_parameters_cellular_potts_class & model,
			 const type_system_class & type_system,
			 const long int & cell_index
			 ) const;
public: void show_cells(
			const model_parameters_cellular_potts_class & model
			) const;
  //
public: void set_buffer_cell(
			     const model_parameters_cellular_potts_class & model,
			     const type_system_class & type_system
			     );
  //
public: int get_type(const long int & cell_index)const;
  //
public: void set_type(
		      const long int & cell_index,
		      const int & type_index
		      );
  //
public: long int get_buffer_cell() const;
  //
public: void calculate_cell_volumes(
				    const model_parameters_cellular_potts_class & model,
				    const std::vector<long int> & configuration,
				    std::vector<long long int> & cell_volumes
				    ) const;
  //
public: void assign_cell_origins(
				 const model_parameters_cellular_potts_class & model,
				 const std::vector<long int> & configuration,
				 std::vector<long long int> & cell_origins
				 ) const;
  //
public: void calculate_cell_total_positions(
					    const model_parameters_cellular_potts_class & model,
					    site_system_class & site_system,
					    const std::vector<long int> & configuration,
					    const std::vector<long long int> & cell_origins,
					    std::vector<long long int> & cell_total_positions
					    ) const;
  //
  //public: const std::vector<long long int> calculate_cell_total_positions(
  //									const model_parameters_cellular_potts_class & model,
  //									site_system_class & site_system,
  //									const std::vector<long int> & configuration,
  //									const std::vector<long long int> & cell_origins
  //									) const;
  //
public: std::vector<double> calculate_cell_position(
						    const model_parameters_cellular_potts_class & model,
						    const site_system_class & site_system,
						    const std::vector<long long int> & cell_origins,
						    const std::vector<long long int> & cell_total_positions,
						    const std::vector<long long int> & cell_volumes,
						    const long int & cell_index
						    ) const;
  //
public: void calculate_cell_position_from_cell_origin(
						      const model_parameters_cellular_potts_class & model,
						      const site_system_class & site_system,
						      const std::vector<long long int> & cell_origins,
						      const std::vector<long long int> & cell_total_positions,
						      const std::vector<long long int> & cell_volumes,
						      const long int & cell_index,
						      std::vector<double> & relative_cell_position 
						      ) const;
  //
public: void calculate_cell_position_in_system(
					       const model_parameters_cellular_potts_class & model,
					       const site_system_class & site_system,
					       const std::vector<long long int> & cell_origins,
					       const std::vector<long long int> & cell_total_positions,
					       const std::vector<long long int> & cell_volumes,
					       const long int & cell_index,
					       std::vector<double> & relative_cell_position 
					       ) const;
  //
public:void  show_cell_positions(
				 const model_parameters_cellular_potts_class & model,
				 const site_system_class & site_system,
				 const std::vector<long int> & configuration,
				 const std::vector<long long int> & cell_origins,
				 const std::vector<long long int> & cell_total_positions,
				 const std::vector<long long int> & cell_volumes
				 )const;
//
public:void check_cell_weight_polarity(
				 const model_parameters_cellular_potts_class & model,
				 site_system_class & site_system,
				 const std::vector<long int> & configuration,
				 const std::vector<long long int> & cell_origins,
				 const std::vector<long long int> & cell_total_positions,
				 const std::vector<long long int> & cell_volumes
				)const;
  //
public:void set_mobile_table();
  //
  /*======================
    Constructor
   =======================*/
public:  cell_system_class(
			   const model_parameters_cellular_potts_class & model,
			   const type_system_class & type_system
			   );
  //
};
#endif // __CELL__

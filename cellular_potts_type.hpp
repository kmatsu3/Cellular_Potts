#ifndef __CELL_TYPE__
#define __CELL_TYPE__
#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
#include "toolbox.hpp"
#include <vector>
//
class type_cellular_potts_class {
  /*======================
    Members
   =======================*/
private: long int Number_of_cells;
private: std::string Species;
private: long int Natural_volume;
private: double Balk_modulus;
private: long int Natural_perimeter;
private: double Persistent_time;
private: double Adhesion_sensitivity;
private: double Polarity_sensitivity;
private: double Quadrapolarity_sensitivity;
private: double Field_sensitivity;
private: std::vector<double> Natural_polarity;
private: std::vector<double> Polar_coupling_table;
private: std::vector<double> Isotropic_Adhesion_Coupling_constant;
private: std::vector<double> Dipolar_Adhesion_Coupling_constant;
private: double Dipolar_Adhesion_Coupling_basal;
private: std::vector<double> Quadrapolar_Adhesion_Coupling_constant;
private: double Quadrapolar_Adhesion_Coupling_basal;
  /*======================
    Methods
   =======================*/
public:std::string set_type_double_value(
					 const std::string & value_name,
					 const double & value
					 );
public:std::string set_type_longint_value(
					  const std::string & value_name,
					  const int & value
					  );
public:std::string set_type_string_value(
					 const std::string & value_name,
					 const std::string & value
					 );
public:std::string add_type_double_vector_value(
						const std::string & value_name,
						const std::vector<double> & value
						);
public: void show() const;
  //
public: long int get_number_of_cells() const;
  /*
public: const std::vector<double> get_adhesion_couplings(const std::string adhesion_type) const;
  */
public: void get_adhesion_couplings(
					  const std::string & adhesion_type,
					  std::vector<double> & coupling_constants
					  ) const;
  //
public: void get_adhesion_basal(
				      const std::string & adhesion_type,
				      double & basal_value
				      ) const;
  //
public: long int get_natural_volume() const;
public: double get_balk_modulus() const;
public: double get_persistent_time() const;
public: double get_adhesion_sensitivity() const;
public: double get_polarity_sensitivity() const;
public: double get_quadrapolarity_sensitivity() const;
public: double get_field_sensitivity() const;
  /*======================
    Constructor
   =======================*/
public:
  type_cellular_potts_class();
  //  explicit type_cellular_potts_class(const type_cellular_potts_class & cell_type);
};
/*==========================
 Wrapping class for type
============================*/
class type_system_class 
{
  /*======================
    Members
   =======================*/
private: std::vector<type_cellular_potts_class> cell_type;
private: int buffer_type;
private: long int long_int_error_number;
private: int int_error_number;
  //private: long int initializer=-1;
  /*======================
    Temporary Members for methods
   =======================*/
private: std::string structure_item;
  /*======================
    Methods
   =======================*/
public:
  void initialize_type(
		       const model_parameters_cellular_potts_class & model
		       );
public: void show_typelist();
public: void show_type(const int & type_index);
public: int get_buffer_type() const;		 
private:void set_structure_type_item(
				     int type_index,
				     std::string child
				     );
  //
public: long int get_number_of_cells(
				     const int & type_index
				     ) const;
  //
  /*
public: const std::vector<double> get_adhesion_couplings(
							 const int type_index, 
							 const std::string adhesion_type
							 ) const;
  */
public: void get_adhesion_couplings(
				    const int & type_index, 
				    const std::string & adhesion_type,
				    std::vector<double> & adhesion_couplings
				    ) const;
  //
public: void get_adhesion_basal(
				const int & type_index, 
				const std::string & adhesion_type,
				double & adhesion_basal
				) const;
  //
public: long int get_natural_volume(const int & type_index) const;
public: std::vector<long int> get_natural_volumes() const;
  //
public: double get_balk_modulus(const int & type_index) const;
public: std::vector<double> get_balk_moduluses() const;
  //
public: double get_persistent_time(const int & type_index) const;
public: void get_persistent_times(
					std::vector<double> & persistence_times
					) const; 
public: void get_adhesion_sensitivities(
					std::vector<double> & adhesion_sensitivities
				      ) const; 
  //
public: double get_adhesion_sensitivity(const int & type_index) const;
  //
public: double get_polarity_sensitivity(const int & type_index) const;
  //
public: double get_quadrapolarity_sensitivity(const int & type_index) const;
  //
public: double get_field_sensitivity(const int & type_index) const;
  //
public: double calculate_cell_scale() const;
  //
public: void set_parameter_double(
				  const int & type_index,
				  const std::string & parameter_identifier,
				  const double & value
				  );
  //
public: void set_parameter_double_vector(
					 const int & type_index,
					 const std::string & parameter_identifier,
					 const std::vector<double> & values
					 );
  //
  /*======================
    Constructor
   =======================*/
public:type_system_class();
};
#endif //__CELL_TYPE__

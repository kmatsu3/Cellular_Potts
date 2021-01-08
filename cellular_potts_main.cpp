#include "cellular_potts_definition.hpp"
#include "cellular_potts_type.hpp"
#include "cellular_potts_cell.hpp"
#include "cellular_potts_site.hpp"
#include "cellular_potts_random.hpp"
#include "cellular_potts_simulation.hpp"
#include "cellular_potts_schedule.hpp"
#include "cellular_potts_observation.hpp"
#include "cellular_potts_shape.hpp"
#include "cellular_potts_region.hpp"
#include "cellular_potts_state.hpp"
#include "cellular_potts_motion_trace.hpp"
int main()
{
  fprintf(stderr,"program is started.\n");
  model_parameters_cellular_potts_class model;
  type_system_class cell_type_system;
  model.input_model_parameters();
  model.show_model_parameters();
  cell_type_system.initialize_type(model);
  cell_type_system.show_typelist();
  adhesion_system_class adhesion_system(model, cell_type_system);
  adhesion_system.initialize_adhesion();;
  adhesion_system.show_adhesion();
  cell_system_class cell_system(model,cell_type_system);
  cell_system.show_cells(model);
  site_system_class site_system(model);
  site_system.generate_sites(model);
  region_system_class region_system(model);
  //
  //  cell_displacement_class cell_track(model);
  cell_displacement_system cell_track(model,cell_type_system);
  //  site_system.show_site_list(model); // for debug
  random_system_seed_class random_seed(model);
  fprintf(stderr,"seed\n");
  random_system_setting_class random_setting(
					     model,
					     random_seed.get_seeds("trial"),
					     random_seed.get_seeds("site_choice"), 
					     random_seed.get_seeds("neighbor_choice"),
					     random_seed.get_seeds("initializer"),
					     model.get_number_of_sites(),
					     site_system.get_number_of_neighbor_sites()
					     );
  fprintf(stderr,"seed set\n");
  simulation_system_class simulation(
				     model, 
				     random_setting.get_engine("trial"),
				     random_setting.get_distribution_double("trial"),
				     random_setting.get_engine("site_choice"),
				     random_setting.get_distribution_longlongint("site_choice"),
				     random_setting.get_engine("neighbor_choice"),
				     random_setting.get_distribution_int("neighbor_choice"),
				     random_setting.get_engine("initializer"),
				     random_setting.get_distribution_double("initializer")
				     );
  fprintf(stderr,"Starting simulation generator\n");
  schedule_system_class schedule(model,cell_type_system,simulation);
  fprintf(stderr,"Starting observation generator\n");
  observation_system_class macro_data(model,cell_system,cell_type_system,simulation);
  observables_type_system_class observables;
  shape_system_class shape_system;
  fprintf(stderr,"Starting state generator\n");
  
  state_system_class state(model,site_system,cell_system,cell_type_system);
  fprintf(stderr,"model constract\n");
  //
  schedule.monte_carlo(
		       model,
		       cell_system,
		       cell_type_system,
		       site_system,
		       region_system,
		       simulation,
		       state,
		       cell_track,
		       macro_data,
		       observables,
		       shape_system,
		       adhesion_system
		       );
  //
  fprintf(stderr,"state updated\n");
  state.finalize_state(model,cell_system,cell_type_system,site_system);
  fprintf(stderr,"program is finished.\n");
  //
}

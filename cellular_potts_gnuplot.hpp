#ifndef IO_CELLULAR_GNUPLOT_
#define IO_CELLULAR_GNUPLOT_
#include <string>
#include <vector>
#include <cstdio>
#include "cellular_potts_io.hpp"
class gnuplot_driver
{
  /*=====================================
    Members for gnuplot format
   =====================================*/
private: std::string plot;
private: std::string xrange;
private: std::string yrange;
private: std::string zrange;
private: std::string xtics;
private: std::string ytics;
private: std::string ztics;
private: std::string xlabel;
private: std::string ylabel;
private: std::string zlabel;
private: std::string size;
private: std::string ratio;
private: std::string view;
private: std::string pm3d;
private: std::string unset_surface;
private: std::string unset_arrow;
private: std::string load;
private: std::string output;
private: std::string file_header;
private: std::string terminal;
private: std::string pause;
  /*=====================================
    Members for terminal format (defined in constractor)
   =====================================*/
private: std::string eps;
private: std::string landscape;
private: std::string png;
private: std::string x_terminal;
private: std::string legend;
private: std::string key;
  /*=====================================
    Members for initializetion
   =====================================*/
private: const int type_max_init;
private: const int type_min_init;
  /*=====================================
    Functions for set
   =====================================*/
public: void set_axis(
		      int direction_index,
		      double min,
		      double max,
		      double tics,
		      std::string label
		      );
  //
public: void set_form(
		      double ratio,
		      double size,
		      std::string view_text,
		      std::string pm3d_text
		      );
public: void make_plot_file(
			    const int & iteration,
			    const std::string & terminal_identifier,
			    const std::string & plot_command,
			    const long int & number_of_steps,
			    const int & number_of_sweeps
			    ) ;
  //
public: std::string get_terminal(
				 const std::string & terminal_identifier
				 ) const;
  //
public: void output_pm3d_data(
			      const std::vector<std::vector<long int> > & data,
			      const long int & number_of_planes,
			      const long int & number_of_steps,
			      const int & number_of_sweeps,
			      const long int & x_dimension,
			      const long int & y_dimension
			      ) const;
//
public: void output_vectors_on_loading_file(
					    const std::vector<std::vector <double> > & coordinates,
					    const std::vector<std::vector <double> > & arrows,
					    const long int & plane_index,
					    const long int & number_of_steps,
					    const int & number_of_sweeps,
					    const int & space_dimension,
					    const double & arrow_size
					    );
  //
public: void output_vectors_with_type_on_loading_file(
						      const std::vector<std::vector <double> > & coordinates,
						      const std::vector<std::vector <double> > & arrows,
						      const std::vector<int> & types,
						      const long int & plane_index,
						      const long int & number_of_steps,
						      const int & number_of_sweeps,
						      const int & space_dimension,
						      const double & arrow_size
						      );
  //
public: gnuplot_driver(std::string input_file_header);
};
#endif // #ifndef IO_CELLULAR_GNUPLOT_

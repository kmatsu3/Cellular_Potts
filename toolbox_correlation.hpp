/ Header for correlation class
#ifndef __CORRELATION__TOOLBOX__
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <cmath>
#include <list>
//
class correlation_evaluation{
  /*=================================
    Definitions of class
   =================================*/
public class element{
protected:std::vector<double> position;
protected:std::vector<double> varaible;
}
public class correlation_system{
public:std::string boundary_condition;
public:int number_of_elements;
public:std::vector<double> size;
public:int number_of_bins;
public:double bin_size;
  /*=================================
    Operators
   =================================*/
  /*=================================
    Constructor & Destractor
   =================================*/
public: void correlation_system();
}
private:long long int number_of_elements;
private:std::vector<double> system_dimension;
private:int number_of_bins;
private:double bin_size;
  /*=================================
    Definitions of general paramerers
   =================================*/
private:vector<int> count;
private:vector<double> relative_vector;
private:vector<double> work_vector;
  /*=================================
    Operators
   =================================*/
public:void calculation(
			const correlation_system & system,
			const std::vector<element> & elements,
			std::vector<double> & correlation
			);
  //
private:void shifted_position(
			      const std::vector<double> & input,
			      std::vector<double> & output,
			      const std::vector<double> & dimensions
			      ) const ;
  //
private:void absolute(
		      const std::vector<double> & input,
		      double & output
		      );
  //
private:void distribute(
			const correlation_system & system,
			const double & relative_length,
			int & bin
			);
  //
private:void init(
		  const correlation_system & system
		  );
  /*=================================
    Constructor & Destractor
   =================================*/
public: void correlation_evaluation(correlation_system);
}
#define __CORRELATION__TOOLBOX__
#endif // __CORRELATION__TOOLBOX__

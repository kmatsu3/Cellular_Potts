// Source file for correlation class
#include "toolbox_correlation.hpp"
void correlation_evaluation::correlation_system::correlation_system{
  boundary_condition = "periodic";
  number_of_element = 0;
};
//
  /*=================================
    Operators
   =================================*/
void correlation_evaluation::calculation(
					 const correlation_system & system,
					 const std::vector<element> & elements,
					 std::vector<double> & correlation
					 )
{
  int number_of_elements = (int)elements.size();
  int space_dimension = (int) elements.position.size();
  double absolute_value;
  int number_of_bins = correlation.size();
  int bin_index;
  for(bin_index=0;bin_index<number_of_bins;bin_index++)
    {
      correlation[bin_index] = 0.0;
    };
  for(int element_index=0;element_index<number_of_elements;element_index++)
    {
      for(int sub_index=element_index+1;sub_index<number_of_elements;sub_index++)
	{
	  for(int space_index=0;space_index<space_dimension;space_index++)
	    {
	      work_vector[space_index]
		=elements[element_index].position[space_index]
		-elements[sub_index].position[space_index];
	    }
	  //
	  shifted_position(
			   work_vector,
			   relative_vector,
			   system.size
			   );
	  //
	  absolute(
		   relative_vector,
		   absolute_value
		   );
	  //
	  distribute(
		     system,
		     absolute_value,
		     bin_index
		     );
	  //
	  if(bin_index<number_of_bins)
	    {
	      count[bin_index]++;
	      for(int space_index=0;space_index<space_dimension;space_index++)
		{
		  correlation[bin_index]
		    +=elements[element_index].variable[space_index]
		    *elements[sub_index].variable[space_index];
		};
	    };
	};
    };
  //
  for(bin_index=0;bin_index<number_of_bins;bin_index++)
    {
      if(count[bin_index]>0)
	{
	  correlation[bin_index]
	    =correlation[bin_index]/(double)count[bin_index];
	};
    };
};
//
void correlation_evaluation::shifted_position(
					      const std::vector<double> & input,
					      std::vector<double> & output,
					      const std::vector<double> & dimensions
					      )
const {
  std::vector<double>::const_iterator input_index=input.begin();
  std::vector<double>::iterator output_index=output.begin();
  std::vector<double>::const_iterator space_index=dimensions.begin();
  int output_index=0;
  while(input_index != input.end())
    {
      if(fabs(*input_index)>fabs(*space_index)) {
	if((*input_index)<0.0)
	  {
	    (*output_index)=(*input_index)+fabs(*space_index);
	  } else {
	    (*output_index)=(*input_index)-fabs(*space_index);
	  
	  };
      };
      input_index++;
      output_index++;
      space_index++;
    };
};
//
void correlation_evaluation::absolute(
				      const std::vector<double> & input,
				      double & output
				      )
{
  output = std::inner_product(input.begin(),input.end(),input.begin(),0.0);
};
//
void correlation_evaluation::distribute(
					const correlation_system & system,
					const double & relative_length,
					int & bin
					);
{
  bin = (int)(relative_length/system.bin_size);
};
//
void correlation_evaluation::init(
				  const std::vector<element> & elements,
				  const std::vector<double> & system_dimensions
				  )
{
  std::vector<element>::const_iterator element_index=elements.begin();
  for
  count.clear();
  for(int bin_index=0;bin_index<system.number_of_bins;bin_index++)
    {
      count.push_back(0.0);
    };
  relative_vector.clear();
  work_vector.clear();
  for(int index=0;index<system.size.size();index++)
    {
      relative_vector.push_back(0.0);
      work_vector.push_back(0.0);
    };
};

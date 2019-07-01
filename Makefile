System_source  = cellular_potts_io.cpp     \
                 cellular_potts_gnuplot.cpp\
                 toolbox.cpp
Modules_source = cellular_potts_definition.cpp     \
                 cellular_potts_random.cpp         \
                 cellular_potts_simulation.cpp     \
                 cellular_potts_polarity.cpp       \
                 cellular_potts_adhesion.cpp       \
                 cellular_potts_type.cpp           \
                 cellular_potts_cell.cpp           \
                 cellular_potts_site.cpp           \
                 cellular_potts_hamiltonian.cpp    \
                 cellular_potts_adhesive_hamiltonian.cpp    \
                 cellular_potts_polarity_motion.cpp\
                 cellular_potts_schedule.cpp  	   \
                 cellular_potts_observation.cpp    \
                 cellular_potts_observation_type.cpp    \
                 cellular_potts_region.cpp         \
                 cellular_potts_motion_trace.cpp   \
                 cellular_potts_state.cpp
Main_source    = cellular_potts_main.cpp
Source = $(System_source)         \
         $(Modules_source)        \
         $(Main_source)
Objects = $(Source:%.cpp=%.o)
System_headers = $(System_source:%.cpp=%.hpp) \
Headers = $(System_source:%.cpp=%.hpp)        \
          $(Modules_source:%.cpp=%.hpp)
Dependences=$(System_source:%.cpp=%.d)  \
	    $(Modules_source:%.cpp=%.d)
Program=cpm
Cc = /opt/intel/bin/icc
Debugs= -g -debug -traceback
Options= -c -static -MMD -MP
ifeq ($(CC),icc_fast)
Cc = /opt/intel/bin/icc
Debugs= 
Options=  -O3 -static -ipo -no-prec-div -c -MMD -MP
endif
ifeq ($(CC),icc_gprof)
Cc = /opt/intel/bin/icc
Debugs= 
Options= -g  -pg -c -MMD -MP
endif
ifeq ($(CC),icc_profgen)
Cc = /opt/intel/bin/icc
Debugs= -O3 -xHOST -O3 -ipo -no-prec-div -prof_gen
Options= -c -MMD -MP
endif
ifeq ($(CC),icc_profuse)
Cc = /opt/intel/bin/icc
Debugs= -O3 -prof_use
Options= -c -static -MMD -MP
endif
ifeq ($(CC),gcc)
Cc = g++
Debugs= -g -ggdb -Wall -pg
Options= -c -static -MMD -MP
endif
ifeq ($(CC),gcc_prof)
Cc = g++
Debugs= -pg
Options= -O0 -c -static -MMD -MP
endif
ifeq ($(CC),gcc_prof5)
Cc = g++
Debugs= -pg
Options= -O5 -c -static -MMD -MP
endif
ifeq ($(CC),gcc_fast)
Cc = g++
Debugs=
Options= -O5 -c -static -MMD -MP
endif
ifeq ($(CC),clang)
Cc = gcc 
LIBs= -I/opt/local/include/ -I/usr/local/include/
Debugs=  -Wall  
Options= -c -static -MMD -MP
endif
-include $(Dependences)
.cpp.o:
	$(Cc) $(Options) $(Debugs) $< -o $@ 
cpm: $(Objects) $(Headers)
	$(Cc) $(Objects) $(Debugs) $(Libs) -v -o $@
cellular_potts_main.o:cellular_potts_main.cpp
cellular_potts_state.o:cellular_potts_state.cpp
cellular_potts_simulation.o:cellular_potts_simulation.cpp
cellular_potts_site.o:cellular_potts_site.cpp
cellular_potts_cell.o:cellular_potts_cell.cpp
cellular_potts_type.o:cellular_potts_type.cpp
cellular_potts_adhesion.o:cellular_potts_adhesion.cpp
cellular_potts_polarity.o:cellular_potts_polarity.cpp
cellular_potts_hamiltonian.o:cellular_potts_hamiltonian.cpp
cellular_potts_adhesive_hamiltonian.o:cellular_potts_adhesive_hamiltonian.cpp
cellular_potts_polarity_motion.o:cellular_potts_polarity_motion.cpp
cellular_potts_schedule.o:cellular_potts_schedule.cpp
cellular_potts_observation.o:cellular_potts_observation.cpp
cellular_potts_observation_type.o:cellular_potts_observation_type.cpp
cellular_potts_random.o:cellular_potts_random.cpp
cellular_potts_region.o:cellular_potts_region.cpp
cellular_potts_region.o:cellular_potts_motion_trace.cpp
cellular_potts_definition.o:cellular_potts_definition.cpp
cellular_potts_gnuplot.o:cellular_potts_gnuplot.cpp
cellular_potts_io.o:cellular_potts_io.cpp
toolbox.o:toolbox.cpp
clean:
	rm -f $(Dependences) $(Objects) $(Program) *~

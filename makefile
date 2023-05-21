
# Path to EZDFFTPACK library.
EZDFFTPTH = ./ezdfftpack

# Command-Line System Options: Linux, Windows.
CMD = Linux
# FC options: mpifort, gfortran.
FC = mpifort

# Executable name.
TRGT = SH23

# .f90 files.
SRCS = $(wildcard *.f90)
# .o files.
OBJS = $(patsubst %.f90, %.o, $(SRCS))
# .mod files.
MODS = $(patsubst %.f90, %.mod, $(SRCS))

# Compiler flags.
FCFLAGS = -O3 -Wall
# Include file search path, e.g., C:\lib\mpi\include.
INCPTHS = $(EZDFFTPTH) 
# Library search path, e.g., C:\lib\mpi\lib.
INCLIBS = $(EZDFFTPTH) 
# Link with specific libraries (.a or .so files), e.g. EZ_PARALLEL.
LIBNAMES = EZDFFTPACK 

# Add appropriate prefixes to include file search paths, library search paths,
# and include libraries.
INCFLAGS = $(addprefix -I,$(INCPTHS))
INCLIBFLAGS = $(addprefix -L,$(INCLIBS))
LKLIBFLAGS = $(addprefix -l,$(LIBNAMES))

# Links objects to the executable.
$(TRGT): $(OBJS)                # Make executable.
	$(FC) $(OBJS) -o $(TRGT) $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

# Make main.o.
main.o: main.f90 iterations.o initialize.o parallel_structs.o parallel_tools.o write_outputs.o
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

# Compile iterations.o
iterations.o: iterations.f90 initialize.o parallel_structs.o parallel_tools.o write_outputs.o
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)	
	
# Complie write_outputs.o
write_outputs.o: write_outputs.f90 initialize.o parallel_structs.o parallel_tools.o
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)	

# Compile initialize.o
initialize.o: initialize.f90 parallel_structs.o parallel_tools.o 
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)	

# Compile all independent .f90 files. 
%.o: %.f90
	$(FC) -c $< $(INCFLAGS) $(INCLIBFLAGS) $(LKLIBFLAGS) $(FCFLAGS)

.PHONY : clean tidy
# Clean deletes everything besides source code.
clean:
	rm -rf *.o *.mod *.csv $(TRGT) 
	rm -rf output_parallel output_serial
	rm -rf ./output_jpg/*.jpg
	rm -rf ./output_jpg/*.eps
	rm -rf ./output_gif/*.gif
	rm -rf ./output_fig/*.fig
	rm -rf ./output_fig/*.jpg
# Tidy deletes all .o and .mod files, keeping the executable.
tidy:
	rm -rf *.o *.mod $(TRGT)






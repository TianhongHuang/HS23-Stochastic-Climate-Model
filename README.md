# Stochastic-Climate-Model
This repo shares the code accomplished simulation of an idealized climate model with spatiotemporal stochastic clouds and planetary-scale circulation, which illustrates study of climate-change problems, changes in cloud cluster statistics, Walker circulations and extreme rainfall events.

# Numerical Setup
1.2: the numerical method used for simulating model SH23 is the operator splitting method.
1.3: the computer code is written in the Fortran programming language with being parallelized using message passing interface (MPI).
1.4: the parallel package used here is build on `EZ_PARALLEL` [distribution](https://github.com/jasonlturner/EZ_PARALLEL_project)

# Before Executing
2.1: all parameters are written in `NAMELIST`, if you want to modify some specific parameters you can go straightforward into that file.
2.2: the eigen pairs data are generated [offline] via Matlab code `eigen_data.m`, before compling executable program you should generate eigen pairs as .dat file in path [./eigen_info], making sure the parameters in `NAMELIST` are consistent with those in `eigen_data.m` (for example number of processors, domain size, etc)

# Running Simulations 
3.1: using command [make] via `makefile` to complie executable program [SH23], ensure that you get the right path to `EZDFFTPACK` library.
3.2: running [SH23] by OpenMPI command [mpirun -np xx ./SH23], '-np xx' is the number of processors we used, for example if we want to use 20 cores to execute program we just need to type `mpirun -np 20 ./SH23`.
3.3: side notes about all .f90 files
    - `parallel_structs.f90`, `parallel_tools.f90`: setting parallel frameworks used by discrete Fourier transform and spectral method.
    - `main.f90`: setting the complie framework of executable program.
    - `initialize.f90`: initializing all parameters and variables used during simulation.
    - `iterations.f90`: describing modules that handling model simulation, involving 4 schemes.
    - `write_outputs.f90`: outputting correspond data as .csv file in [./output_parallel].

# Plotting Figures
4.1: all outputs data are stored as .csv file in path [./output_parallel], categoried via processor id and variable names.

## Collecting Data
4.2: using Matlab code to load all .csv file and combine data from each individual processors, running `READING_xx.m` in the same folder of executable program, it would store useful data as .mat file for further plotting work.
    - `READING_varbs.m`: reading and storing time-averaged and time-meridionally-averaged data as [stats_varbs_ast.mat].
    - `READING_hist_rain.m`: reading and storing histogram data of rainfall events as [hist_rain.mat].
    - `READING_hist_cluster.m`: reading and storing histogram data of clusters' size as [hist_cluster.mat].

## Visualizing Statistics
4.3: adding right tags to all .mat files, for example, 
    - `stats_varbs_11.mat`: model variables from standard simulation.
    - `stats_varbs_12.mat`: model variables from climate-change simulation.
    - `stats_varbs_1.mat` to `stats_varbs_5.mat`: model variables from simulations with different cloud albedo values, used in sensitivity studies.
    - using the similar fashion for adding tags to `hist_rain.mat` and `hist_cluster.mat` too.
4.4: puting all .mat files with right tags together with `FIG_xxx.m` file to make comparison figures,
    - `FIG_climate_change.m`: figure that comparing time-meridionally-averaged model variables between standard simulation and climate-change simulation.
    - `FIG_histogram.m`: figure that comparing histogram data between standard simulation and climate-change simulation.
    - `FIG_sens_study.m`: figure that comparing time-meridionally-averaged model variables and histogram data from simulations with different cloud albedo values (Ab or Af).
    
# Contact Information
If you have any comments, questions, concerns, or ideas about this program, please contact Tianhong Huang at thuang76@wisc.edu
    

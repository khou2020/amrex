## Parallel I/O Kernel Case Study For AMRex

This benchmark program is a case study of parallel I/O kernel from the
[AMReX](https://github.com/AMReX-Codes/amrex) software framework library. The
benchmark program, can be used to evaluate
[PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF) library for its
performance on the I/O patterns used by AMReX plot file.

This benchmark follows the flow of the HDF5 benchmark. It wrties the same plot file data as the HDF5 benchmark except 
that the file are in NetCDF instead of HDF5. Aside from the NetCDF plot file. The benchmark also output the plot file 
in AMRex's original file format. It also write checkpoint files in raw AMRex format. Those non-NetCDF files are not measured 
by the benchmark program. We keep it to maintain consistancy to the HDF5 benchmark.

The AMRex grid is a 3-dimensional array of size N * N * N where N is configurable via option. 
It is divided into rectangular subdomain which is distributed among processes.
To generate the plot file, the processes gather data from subdomain it is responsible for into a continuous buffer. 
The data form each process is then stacked together to be stored in a single 1 dimensional variable by calling 
ncmpi_put_vara_double either collectively or independently based on compile option.

The benchmark program writes 4 NetCDF variables. The first one is a 1-dimensional integer variable of size M recording 
the rank of process responsible to individual subdomain. The second one is a 2 dimensional variable of size m * 6 
representing the bounding box of each subdomain. Each row represents 1 subdomain. The third variable is a 1-dimensional 
64 bit integer variable of size M recording the offset of the subdomain data in the data variable. The fourth variable is 
a one dimensional double variable that store the values in the grid. The size is equal to the size of grid tims number 
of components per cell. Beside variables, the benchmark program also add attributess to describe the file such as grid size.
The benchmark program is implementated using PnetCDF blocking vara APIs. In addition to end to end timings and I/O bandwidths, 
the benchmark also reports time spent in individual steps. 

* Compile command:
  * Edit `GNUmakefile` to customize the compiler, compile options, location of
    PnetCDF library, etc. Variable "PNETCDF_HOME" must be set to the install directory 
    of PnetCDF library.
  * By default, the bnechmark uses collective I/O. To use independent I/O, add `-NCINDEP` 
    to DEFINES in the `GNUmakefile`.
  * The minimum PnetCDF version tested on the benchmark is 1.10.0.
  * Run command `make` to generate the executable program. The program is named using 
    default AMRex package naming rule. The name of compiled benchmark is `main3d.<options>.ex`
    where <options> varies on different compiliation options.

* Run command:
  * Example run command using `mpiexec` and 16 MPI processes (assuming the compiled program is
    `main3d.gnu.TPROF.MPI`, the configuration file is `inputs`):
    `mpiexec -n 16 ./main3d.gnu.TPROF.MPI.ex inputs`
  * The options are set inside a configuration file that is passed as the first argument to the 
    benchmark program. 
  * Options suppoorted in the configuration file are:
    * ncells
      The size along each dimension of the 3-dimensional grid. The grid will be a 3-dimensional 
      array of size ncells*ncells*ncells.
    * max_grid_size
      Maximum allowable size of each subdomain in the problem domain. The grid is divided into 
      many sub-grid (sub-array) to be distributed among processes. The size of the subarray should 
      not exceed max_grid_size*max_grid_size*max_grid_size.
    * nlevs
      Number of levels of thew grid system. Currently, only 1 level is supported.
    * ncomp
      Number of components in the multifabs. It affects the number of values per grid cell in the plot file.
      The size of the plot file is linearly propotional to ncomp.
    * nppc
      Number of particles in each of the grid cell. It does not affect the size of the plot file.
  * An exaample of the configuration file can be found at ./inputs
    ```
      # Domain size
      ncells = 128 
      # Maximum allowable size of each subdomain in the problem domain; 
      # this is used to decompose the domain for parallel calculations.
      max_grid_size = 32
      # Number of levels
      nlevs = 1
      # Number of components in the multifabs
      ncomp = 6 
      # Number of particles per cell
      nppc = 2 
    ```
  * An example batch script file for running a job with 32 haswell nodes, 32 MPI
    processes per node, on Cori @NERSC is provided in `./slurm.haswell`. 
    It assums the benchmark program is called `amrex_bench` and the configuration file is called `inputs`.

* Example outputs on screen
  ```
    mpiexec -np 4 ./main3d.gnu.TPROF.MPI.ex inputs

    MPI initialized with 4 MPI processes
    AMReX (19.01-75-g807b01a19821) initialized
    Write_NCHEADER_time_1_2_0 = 0.00633597374  seconds.
    Write_NCATT_time_1_3 = 0.0002799034119  seconds.
    Write_NCINIT_time_1 = 0.04796266556  seconds.
    Error at line 245 in WritePlotfilePnetCDF.cpp: Attribute value is inconsistent among processes.
    Error at line 245 in WritePlotfilePnetCDF.cpp: Attribute value is inconsistent among processes.
    Error at line 245 in WritePlotfilePnetCDF.cpp: Attribute value is inconsistent among processes.
    Error at line 245 in WritePlotfilePnetCDF.cpp: Attribute value is inconsistent among processes.
    Write_NCATT_time_2 = 0.0004410743713  seconds.
    Write_NCATT_time_3 = 0.00026512146  seconds.
    Write_NCVAR_time_4 = 0.009718418121  seconds.
    Write_NCATT_time = 0.05839133263  seconds.
    ::---- calling NCPUT for the grid data on level 0
    Write_NCATT_time_5 = 0.01653146744  seconds.
    Write_NCVARPUT_time = 0.1496036053  seconds.
    Write_NCVARPUT_time_since = 1547370031
    Write_NC_time_7_closing = 1.907348633e-06  seconds.
    Write_PNETCDF_time = 0.2281308174  seconds.


    TinyProfiler total time across processes [min...avg...max]: 3.652 ... 3.66 ... 3.667

    ----------------------------------------------------------------------------------------------------
    Name                                                 NCalls  Excl. Min  Excl. Avg  Excl. Max   Max %
    ----------------------------------------------------------------------------------------------------
    ParticleContainer<NSR, NSI, NAR, NAI>::InitRandom()       1      1.275      1.362      1.426  38.88%
    ParticleContainer::RedistributeCPU()                      1     0.7592     0.8138     0.8899  24.27%
    ParticleContainer::RedistributeMPI()                      1     0.3708     0.4228     0.4881  13.31%
    RedistributeMPI_copy                                      1     0.2566     0.2756     0.2847   7.76%
    ParticleContainer::WriteParticles()                       1     0.1609     0.1687     0.1802   4.91%
    RD:convertFromNativeDoubleFormat                         16     0.1474     0.1489     0.1505   4.10%
    NCVarPutGrids                                             1    0.03386     0.1143     0.1496   4.08%
    NCVarPut                                                  1     0.0123    0.04763     0.1282   3.50%
    PD_convert                                               64    0.08452    0.08575    0.08666   2.36%
    ParticleContainer::Checkpoint()                           1    0.05702    0.06147    0.06845   1.87%
    WriteMultiLevelPlotfilePNETCDF                            1    0.02543    0.04502    0.06639   1.81%
    RedistributeMPI_locate                                    1     0.0404    0.04645    0.06028   1.64%
    VisMF::Write(FabArray)                                    1    0.01487    0.01741    0.01852   0.51%
    FabArray::setVal()                                        1     0.0118    0.01495    0.01751   0.48%
    VisMF::WriteHeader                                        1      2e-06    0.00261    0.01044   0.28%
    DistributionMapping::LeastUsedCPUs()                      1    6.8e-05   0.004072   0.008242   0.22%
    WriteMultiLevelPlotfile()                                 1   0.002976   0.004116   0.007523   0.21%
    VisMF::CalculateMinMax                                    1   0.002575   0.002852   0.003164   0.09%
    WriteGenericPlotfileHeader()                              0          0   0.000125     0.0005   0.01%
    FABio_binary::write_header                               48   0.000227  0.0003058    0.00048   0.01%
    FABio::write_header                                      48   0.000113  0.0001495   0.000238   0.01%
    VisMF::FindOffsets                                        1      4e-06    4.6e-05   0.000167   0.00%
    DistributionMapping::SFCProcessorMapDoIt()                1    4.2e-05      6e-05    6.6e-05   0.00%
    VisMF::Header                                             1      2e-06      4e-06      7e-06   0.00%
    PC<NNNN>::Checkpoint:unlink                               0          0   1.25e-06      5e-06   0.00%
    ----------------------------------------------------------------------------------------------------

    ----------------------------------------------------------------------------------------------------
    Name                                                 NCalls  Incl. Min  Incl. Avg  Incl. Max   Max %
    ----------------------------------------------------------------------------------------------------
    ParticleContainer<NSR, NSI, NAR, NAI>::InitRandom()       1      2.901       2.92      2.937  80.10%
    ParticleContainer::RedistributeCPU()                      1      1.509      1.559      1.634  44.55%
    ParticleContainer::RedistributeMPI()                      1      0.671     0.7448      0.833  22.72%
    ParticleContainer::Checkpoint()                           1     0.4615     0.4648     0.4699  12.81%
    ParticleContainer::WriteParticles()                       1     0.3934     0.4034     0.4129  11.26%
    RedistributeMPI_copy                                      1     0.2566     0.2756     0.2847   7.76%
    RD:convertFromNativeDoubleFormat                         16     0.2325     0.2347     0.2372   6.47%
    WriteMultiLevelPlotfilePNETCDF                            1     0.1873      0.207     0.2283   6.23%
    NCVarPut                                                  1     0.1619     0.1619     0.1621   4.42%
    NCVarPutGrids                                             1    0.03386     0.1143     0.1496   4.08%
    PD_convert                                               64    0.08452    0.08575    0.08666   2.36%
    RedistributeMPI_locate                                    1     0.0404    0.04645    0.06028   1.64%
    WriteMultiLevelPlotfile()                                 1    0.02447    0.02762    0.03679   1.00%
    VisMF::Write(FabArray)                                    1    0.02149    0.02338    0.02877   0.78%
    FabArray::setVal()                                        1     0.0118    0.01495    0.01751   0.48%
    VisMF::WriteHeader                                        1      2e-06    0.00261    0.01044   0.28%
    DistributionMapping::SFCProcessorMapDoIt()                1   0.000134   0.004132   0.008308   0.23%
    DistributionMapping::LeastUsedCPUs()                      1    6.8e-05   0.004072   0.008242   0.22%
    VisMF::CalculateMinMax                                    1   0.002575   0.002852   0.003164   0.09%
    FABio_binary::write_header                               48   0.000341  0.0004553   0.000718   0.02%
    VisMF::FindOffsets                                        1      4e-06  0.0001468    0.00057   0.02%
    WriteGenericPlotfileHeader()                              0          0   0.000125     0.0005   0.01%
    FABio::write_header                                      48   0.000113  0.0001495   0.000238   0.01%
    VisMF::Header                                             1      2e-06      4e-06      7e-06   0.00%
    PC<NNNN>::Checkpoint:unlink                               0          0   1.25e-06      5e-06   0.00%
    ----------------------------------------------------------------------------------------------------

    AMReX (19.01-75-g807b01a19821) finalized

  ```
* Output files
  * Each run of the benchmark produces an output netCDF files named
    `plt00000.nc` representing the plot file. If the file already exist, it will be overwritten.
  * The benchmark also output an folder containing the plot file and checkpoint file in raw AMRex format. 
    It is not used in performance evaluation of PnetCDF but was kept to maintain consistancy to the HDF5 benchmark.
  * Here is the netcdf file header of the plot file:
    ```
      ncmpidump plt00000.nc -h
      
      netcdf plt00000 {
      // file format: CDF-5 (big variables)
      dimensions:
              procdataspace_0 = 64 ;
              boxdataspace_0 = 64 ;
              boxdataspace_1 = 6 ;
              offsetdataspace_0 = 65 ;
              dataspace = 4194304 ;
      variables:
              int __level_0_Processors(procdataspace_0) ;
              int __level_0_boxes(boxdataspace_0, boxdataspace_1) ;
              int64 __level_0_data:offsets=0(offsetdataspace_0) ;
              double __level_0_data:datatype=0(dataspace) ;

      // global attributes:
                      :Chombo_globalSpaceDim = 3 ;
                      :Chombo_globaltestReal = 0. ;
                      :Chombo_globaltestString = "vMString::testString" ;
                      :_num_components = 2 ;
                      :_num_levels = 1 ;
                      :_component_0 = "component_0" ;
                      :_component_1 = "component_1" ;
                      :_filetype = "VanillaAMRFileType" ;
                      :_testString = "vMString::testString" ;
                      :__level_0_data_attributescomps = 2 ;
                      :__level_0_dt = 0. ;
                      :__level_0_dx = 0.0078125 ;
                      :__level_0_time = 0. ;
                      :__level_0_prob_domain = 0, 0, 0, 127, 127, 127 ;
                      :__level_0_ref_ratio = 1 ;
      }
    ```
  * Performance evaluation:
    * Configurations:
      1024 processes (32 per haswell node) on Cori @ NERSC
      Lustre 64 stripe with 1 MiB stripe size
      Input file (49 GiB total size):
        ncells = 1024 
        max_grid_size = 32
        nlevs = 1
        ncomp = 6 
        nppc = 2 
    * End to end I/O time (sec):
        Round   1   2   3
        PnetCDF 17  16  28

## Questions/Comments:
email: khl7265@eecs.northwestern.edu

Copyright (C) 2019, Northwestern University.

See [COPYRIGHT](COPYRIGHT) notice in top-level directory.


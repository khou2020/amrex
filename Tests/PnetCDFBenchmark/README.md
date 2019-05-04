## AMRex I/O kernel using PnetCDF

This benchmark program is designed to evaluate [AMReX](https://github.com/AMReX-Codes/amrex) software framework library I/O performance when the I/O module is
implemented using [PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF) library.

This benchmark follows the flow of the HDF5 benchmark. It wrties the same plot file data as the HDF5 benchmark except 
that the file are in NetCDF (https://www.unidata.ucar.edu/software/netcdf/) instead of HDF5. Aside from the NetCDF plot file. 
The benchmark also output the plot file in AMRex's original file format. It also write checkpoint files in raw AMRex format. 
Those non-NetCDF files are not measured by the benchmark program. We keep it to maintain consistancy to the HDF5 benchmark.

* I/O pattern:
  * Data partition:
    * The AMRex grid is a 3-dimensional array of size N * N * N where N is configurable via option. 
    * It is divided into rectangular subdomain which is distributed among processes.
    * The grid is divided into sub-grid, each one is handled by a process
  * Data objects:
    * The benchmark program writes 4 NetCDF variables. 
      * 1-dimensional integer variable of size M recording the rank of process responsible to individual subdomain. 
      * 2 dimensional variable of size m * 6 representing the bounding box of each subdomain. Each row represents 1 subdomain.
      * 1-dimensional 64 bit integer variable of size M recording the offset of the subdomain data in the data variable. 
      * 1-dimensional double variable that store the values in the grid. The size is equal to the size of grid tims number of components per cell. 
    * Beside variables, the benchmark program also add attributess to describe the file such as grid size.

* Implementation:
  * The benchmark program is implementated using PnetCDF blocking vara APIs.
  * In the plot file, the subgrids are flattened and serialized into a 1-dimensional array of bytes.
  * Subgrid data are simply stacked one after another.
  * Processes gather data from subgrid it is responsible for into a continuous buffer.
  * Affter synchronizing the offset anbd size of each subgrid, processes write the data to the variable.
  * The program uses ncmpi_put_vara_double to write the variable data either collectively or independently based on compile option.
  * In addition to end to end timings and I/O bandwidths, the benchmark also reports time spent in individual steps. 

* Building PnetCDF for AMRex benchmark:
  * The benchmark program does not use any advanced feature of PnetCDF, so the default build option will work
  * Download PnetCDF
    git clone https://github.com/Parallel-NetCDF/PnetCDF.git
  * Compile PnetCDF
    ./configure --prefix=<install directory>
    make install
  * For details about building PnetCDF, please refer to https://github.com/Parallel-NetCDF/PnetCDF/blob/master/INSTALL

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
    Write_NCHEADER_time_1_2_0 = 0.009748935699  seconds.
    Write_NCATT_time_1_3 = 0.0004329681396  seconds.
    Write_NCINIT_time_1 = 0.06093502045  seconds.
    Write_NCATT_time_2 = 0.0003917217255  seconds.
    Write_NCATT_time_3 = 0.0007684230804  seconds.
    Write_NCVAR_time_4 = 0.009308576584  seconds.
    Write_NCATT_time = 0.07141041756  seconds.
    ::---- calling NCPUT for the grid data on level 0
    Write_NCATT_time_5 = 0.01284003258  seconds.
    Write_NCVARPUT_time = 0.09048891068  seconds.
    Write_NCVARPUT_time_since = 1547596225
    Write_NC_time_7_closing = 2.145767212e-06  seconds.
    Write_PNETCDF_time = 0.1815581322  seconds.


    TinyProfiler total time across processes [min...avg...max]: 3.735 ... 3.747 ... 3.758

    ----------------------------------------------------------------------------------------------------
    Name                                                 NCalls  Excl. Min  Excl. Avg  Excl. Max   Max %
    ----------------------------------------------------------------------------------------------------
    ParticleContainer<NSR, NSI, NAR, NAI>::InitRandom()       1      1.268      1.286      1.306  34.76%
    ParticleContainer::RedistributeCPU()                      1     0.7699      0.872      0.949  25.25%
    ParticleContainer::RedistributeMPI()                      1     0.3839     0.4785     0.5845  15.55%
    RedistributeMPI_copy                                      1     0.2507     0.2969     0.3209   8.54%
    RD:convertFromNativeDoubleFormat                         16     0.1662     0.1757     0.1923   5.12%
    ParticleContainer::WriteParticles()                       1     0.1731     0.1794     0.1893   5.04%
    PD_convert                                               64    0.09295    0.09433    0.09741   2.59%
    NCVarPutGrids                                             1    0.04325    0.07374    0.09048   2.41%
    WriteMultiLevelPlotfilePNETCDF                            1    0.03166     0.0519    0.07936   2.11%
    ParticleContainer::Checkpoint()                           1    0.07165    0.07267    0.07439   1.98%
    RedistributeMPI_locate                                    1    0.05693    0.06287    0.06847   1.82%
    NCVarPut                                                  1    0.01209    0.02886    0.05941   1.58%
    VisMF::Write(FabArray)                                    1    0.01751    0.01987    0.02081   0.55%
    FabArray::setVal()                                        1    0.01203    0.01493    0.01788   0.48%
    WriteMultiLevelPlotfile()                                 1   0.003566   0.006393    0.01486   0.40%
    VisMF::WriteHeader                                        1      2e-06    0.00327    0.01308   0.35%
    VisMF::CalculateMinMax                                    1   0.002674    0.00329   0.003667   0.10%
    DistributionMapping::LeastUsedCPUs()                      1   0.000864  0.0008698   0.000876   0.02%
    FABio_binary::write_header                               48   0.000205   0.000417   0.000661   0.02%
    WriteGenericPlotfileHeader()                              0          0  0.0001555   0.000622   0.02%
    VisMF::FindOffsets                                        1      5e-06   0.000149   0.000577   0.02%
    DistributionMapping::SFCProcessorMapDoIt()                1   0.000301  0.0003447    0.00036   0.01%
    FABio::write_header                                      48    9.9e-05  0.0001892   0.000359   0.01%
    PC<NNNN>::Checkpoint:unlink                               0          0   1.75e-06      7e-06   0.00%
    VisMF::Header                                             1      3e-06   3.75e-06      4e-06   0.00%
    ----------------------------------------------------------------------------------------------------

    ----------------------------------------------------------------------------------------------------
    Name                                                 NCalls  Incl. Min  Incl. Avg  Incl. Max   Max %
    ----------------------------------------------------------------------------------------------------
    ParticleContainer<NSR, NSI, NAR, NAI>::InitRandom()       1      2.971      2.996      3.019  80.34%
    ParticleContainer::RedistributeCPU()                      1      1.679       1.71      1.742  46.34%
    ParticleContainer::RedistributeMPI()                      1     0.7397     0.8383     0.9718  25.86%
    ParticleContainer::Checkpoint()                           1     0.5157     0.5221     0.5308  14.12%
    ParticleContainer::WriteParticles()                       1     0.4426     0.4494     0.4591  12.22%
    RedistributeMPI_copy                                      1     0.2507     0.2969     0.3209   8.54%
    RD:convertFromNativeDoubleFormat                         16     0.2636       0.27      0.286   7.61%
    WriteMultiLevelPlotfilePNETCDF                            1     0.1342     0.1545     0.1819   4.84%
    NCVarPut                                                  1     0.1026     0.1026     0.1027   2.73%
    PD_convert                                               64    0.09295    0.09433    0.09741   2.59%
    NCVarPutGrids                                             1    0.04325    0.07374    0.09048   2.41%
    RedistributeMPI_locate                                    1    0.05693    0.06287    0.06847   1.82%
    WriteMultiLevelPlotfile()                                 1    0.02785    0.03374    0.05034   1.34%
    VisMF::Write(FabArray)                                    1    0.02429    0.02719    0.03487   0.93%
    FabArray::setVal()                                        1    0.01203    0.01493    0.01788   0.48%
    VisMF::WriteHeader                                        1      2e-06    0.00327    0.01308   0.35%
    VisMF::CalculateMinMax                                    1   0.002674    0.00329   0.003667   0.10%
    DistributionMapping::SFCProcessorMapDoIt()                1   0.001176   0.001215   0.001236   0.03%
    VisMF::FindOffsets                                        1      5e-06  0.0003097    0.00122   0.03%
    FABio_binary::write_header                               48   0.000304  0.0006063    0.00102   0.03%
    DistributionMapping::LeastUsedCPUs()                      1   0.000864  0.0008698   0.000876   0.02%
    WriteGenericPlotfileHeader()                              0          0  0.0001555   0.000622   0.02%
    FABio::write_header                                      48    9.9e-05  0.0001892   0.000359   0.01%
    PC<NNNN>::Checkpoint:unlink                               0          0   1.75e-06      7e-06   0.00%
    VisMF::Header                                             1      3e-06   3.75e-06      4e-06   0.00%
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
              int level_0.Processors(procdataspace_0) ;
              int level_0.boxes(boxdataspace_0, boxdataspace_1) ;
              int64 level_0.data:offsets=0(offsetdataspace_0) ;
              double level_0.data:datatype=0(dataspace) ;

      // global attributes:
                      :Chombo_globalSpaceDim = 3 ;
                      :Chombo_globaltestReal = 0. ;
                      :Chombo_globaltestString = "vMString::testString" ;
                      :num_components = 2 ;
                      :num_levels = 1 ;
                      :component_0 = "component_0" ;
                      :component_1 = "component_1" ;
                      :filetype = "VanillaAMRFileType" ;
                      :testString = "vMString::testString" ;
                      :level_0.data_attributes.comps = 2 ;
                      :level_0.data_attributes.ghost = 0, 0, 0 ;
                      :level_0.dt = 0. ;
                      :level_0.dx = 0.0078125 ;
                      :level_0.time = 0. ;
                      :level_0.prob_domain = 0, 0, 0, 127, 127, 127 ;
                      :level_0.ref_ratio = 1 ;
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


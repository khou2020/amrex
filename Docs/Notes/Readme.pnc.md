# Summary

This AMReX distribution contains a prototype implementation of the PnetCDF I/O module.
The PnetCDF module add methods to access NetCDF formatted file to AMReX data structure class
The module currently supports checkpointing multifab (checkpoint file), checkpointing particle container (particle file), and plot file.

## Enable PnetCDF I/O module
* Enable PnetCDF in GNUmakefile
  + Set "USE_PNETCDF" to "TRUE" in GNUmakefile
* To use the PnetCDF I/O module, you need a PnetCDF library
  + Set "PNETCDF_HOME" to PnetCDF install path in GNUmakefile
* MPI must be enabled to use PnetCDF I/O module
  + Set "USE_MPI" to "TRUE" in GNUmakefile
* PnetCDF does not support OMP
  + Set "USE_OMP" to "FALSE" in GNUmakefile
* The module supports different data layout in NetCDF file
  + 0: Single variable
    + A large variable storing all the grids within the multifab
    + The dimension is AMREX_SPACEDIM + 1. Grids stacks along the first dimension
      + Grid structure preserved
    + Remaining dimensions are made large enough to fit the largest grid along each dimension
      + Can waste space if grid size differs a lot
  + 1: Single variable flattened
    + A large 1-dimensional variable storing all the grids within the multifab
    + Grid data is flattened into a byte stream
      + Structure information lost
    + Grid data appends one after another in the data variable
      + additional metadata variable storing the size of each grid so the location can be calculated
  + 2: Variable per grid
    + One variable per grid
      + Variable size set to exact size of the grid
      + No space waste
      + Data structure preserved
    + Additional metadata variable to link grid variables to multifab
  + Set "PNC_MF_LAYOUT" to the layout id in GNUmakefile
* Example
  ```
  AMREX_HOME ?= ../../
  DEBUG    = TRUE
  DIM    = 3
  COMP    = gnu
  PRECISION = DOUBLE
  TINY_PROFILE = TRUE
  EBASE     = main
  include $(AMREX_HOME)/Tools/GNUMake/Make.defs
  include ./Make.package
  include $(AMREX_HOME)/Src/Base/Make.package
  include $(AMREX_HOME)/Src/Particle/Make.package
  include $(AMREX_HOME)/Tools/GNUMake/Make.rules

  ###################################################
  # Use PnetCDF I/O module
  ################################################### 
  PNETCDF_HOME   = /Path/to/PnetCDF/installation
  USE_PNETCDF = TRUE
  USE_MPI   = TRUE
  USE_OMP   = FALSE
  PNC_MF_LAYOUT=2
  ###################################################
  ```

## Checkpoint file
* VisMF::Write_PNC(ncid, mf, name, setghost)
  + Write a multifab object to an opened NetCDF file
  + User need to open/create the file before calling VisMF::Write_PNC
  + Arguments
    + int ncid: id of the opened NetCDF file. 
    + MultiFab& mf: Reference to the multifab to write
    + string& name: Name of the multifab
    + bool setghost: Whether to fill ghost cell before writing
* VisMF::Write_PNC(fname, mf, name, setghost)
  + An overload of VisMF::Write_PNC(ncid, mf, name, setghost) that takes filename instead of an opened file
  + If the file exists, the mf will be added to the file
  + Arguments
    + string& fname: Name of NetCDF file to write. 
    + MultiFab& mf: Reference to the multifab to write to
    + string& name: Name of the multifab
    + bool setghost: Whether to fill ghost cell before writing
* VisMF::Read_PNC(ncid, mf, name, allow_empty_mf)
  + Read a multifab object from an opened NetCDF file
  + User need to open/create the file before calling VisMF::Read_PNC
  + Arguments
    + int ncid: id of the opened NetCDF file. 
    + MultiFab& mf: Reference to the multifab to store data read
    + string& name: Name of the multifab in the file
    + bool allow_empty_mf: Not used. Added to be consistant to VisMF::Read
* VisMF::Read_PNC(fname, mf, name, allow_empty_mf)
  + An overload of VisMF::Read_PNC(ncid, mf, name, allow_empty_mf) that takes filename instead of an opened file
  + Arguments
    + string& fname: Name of NetCDF file to read from
    + MultiFab& mf: Reference to the multifab to store data read
    + string& name: Name of the multifab in the file
    + bool allow_empty_mf: Not used. Added to be consistant to VisMF::Read
* Example
    ```
    int ncid;
    MultiFab mf;
    string name = "example"

    ncmpi_create(MPI_COMM_WORLD, "checkpoint.nc", 0, MPI_INFO_NULL, &ncid);
    ncmpi__enddef(ncid, 10485760, 0, 0, 0);
    VisMF::Write_PNC(ncid, mf, name, false);
    ncmpi_close(ncid);

    ncmpi_open(MPI_COMM_WORLD, "checkpoint.nc", 0, MPI_INFO_NULL, &ncid);
    VisMF::Read_PNC(ncid, mf, name, false);
    ncmpi_close(ncid);
    ```

## Particle file
* ParticleContainer::Checkpoint_PNC(ncid, name, is_checkpoint, real_comp_names, int_comp_names)
  + Save the ParticleContainer to an opened NetCDF file
  + User need to open/create the file before calling ParticleContainer::Checkpoint_PNC
  + Arguments
    + int ncid: id of the opened NetCDF file. 
    + string& name: Name of the ParticleContainer
    + bool is_checkpoint: If it is for checkpoint
    + string& real_comp_names: Name of real valued parameters
    + string& int_comp_names: Name of integer valued parameters
* ParticleContainer::Restart_PNC(ncid, name)
  + Load the ParticleContainer from an opened NetCDF file
  + User need to open/create the file before calling ParticleContainer::Checkpoint_PNC
  + Arguments
    + int ncid: id of the opened NetCDF file. 
    + string& name: Name of the ParticleContainer
* Example
    ```
    int ncid;
    MyPC myPC(geom, dmap, ba, ref_ratio);
    string name = "example"

    ncmpi_create(MPI_COMM_WORLD, "checkpoint.nc", 0, MPI_INFO_NULL, &ncid);
    ncmpi__enddef(ncid, 10485760, 0, 0, 0);
    myPC.Checkpoint_PNC(ncid, name, true, particle_realnames, particle_intnames);
    ncmpi_close(ncid);

    ncmpi_open(MPI_COMM_WORLD, "checkpoint.nc", 0, MPI_INFO_NULL, &ncid);
    myPC.Restart_PNC(ncid, name);
    ncmpi_close(ncid);
    ```

## Limitation
* Setting communicator is not supported.
  + MPI_COMM_WORLD is assume
    + All process need to participate when using the PnetCDF I/O module

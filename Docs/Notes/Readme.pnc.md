# Summary

A prototype implementation of PnetCDF I/O module is contained in this AMReX distribution
The PnetCDF module add methods to access NetCDF formated file to AMReX datastructure class
The module currently support checkpointing multifab (checkpoint file), checkpointing particle container (particle file), and plot file.

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
  + 0:
  + 1:
  + 2:
  + Set "PNC_MF_LAYOUT" to the layout id in GNUmakefile
* Example
  ```
  AMREX_HOME ?= ../../
  DEBUG	= TRUE
  DIM	= 3
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
  + An overload of VisMF::Write_PNC(ncid, mf, name, setghost) that takes file name instead of opened file
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
  + An overload of VisMF::Read_PNC(ncid, mf, name, allow_empty_mf) that takes file name instead of opened file
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

## Known issues
* All process need to participate when using the PnetCDF I/O module
  + MPI_COMM_WORLD is assume
  + There is no way to set communicator so far

#include <fstream>
#include <iomanip>

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#if defined(BL_USE_MPI) && defined(BL_USE_PNETCDF)
#include <pnetcdf.h>
#endif

namespace amrex {

std::string LevelPath (int level, const std::string &levelPrefix)
{
    return Concatenate(levelPrefix, level, 1);  // e.g., Level_5
}

std::string MultiFabHeaderPath (int level,
                                const std::string &levelPrefix,
                                const std::string &mfPrefix)
{
    return LevelPath(level, levelPrefix) + '/' + mfPrefix;  // e.g., Level_4/Cell
}

std::string LevelFullPath (int level,
                           const std::string &plotfilename,
                           const std::string &levelPrefix)
{
    std::string r(plotfilename);
    if ( ! r.empty() && r.back() != '/') {
	r += '/';
    }
    r += LevelPath(level, levelPrefix);  // e.g., plt00005/Level_5
    return r;
}

std::string MultiFabFileFullPrefix (int level,
                                    const std::string& plotfilename,
                                    const std::string &levelPrefix,
                                    const std::string &mfPrefix)
{
    std::string r(plotfilename);
    if ( ! r.empty() && r.back() != '/') {
	r += '/';
    }
    r += MultiFabHeaderPath(level, levelPrefix, mfPrefix);
    return r;
}


void
PreBuildDirectorHierarchy (const std::string &dirName,
                           const std::string &subDirPrefix,
                           int nSubDirs, bool callBarrier)
{
  UtilCreateCleanDirectory(dirName, false);  // ---- dont call barrier
  for(int i(0); i < nSubDirs; ++i) {
    const std::string &fullpath = LevelFullPath(i, dirName);
    UtilCreateCleanDirectory(fullpath, false);  // ---- dont call barrier
  }

  if(callBarrier) {
    ParallelDescriptor::Barrier();
  }
}


void
WriteGenericPlotfileHeader (std::ostream &HeaderFile,
                            int nlevels,
                            const Vector<BoxArray> &bArray,
                            const Vector<std::string> &varnames,
                            const Vector<Geometry> &geom,
                            Real time,
                            const Vector<int> &level_steps,
                            const Vector<IntVect> &ref_ratio,
                            const std::string &versionName,
                            const std::string &levelPrefix,
                            const std::string &mfPrefix)
{
        BL_PROFILE("WriteGenericPlotfileHeader()");

        BL_ASSERT(nlevels <= bArray.size());
        BL_ASSERT(nlevels <= geom.size());
        BL_ASSERT(nlevels <= ref_ratio.size()+1);
        BL_ASSERT(nlevels <= level_steps.size());

        int finest_level(nlevels - 1);

	HeaderFile.precision(17);

	// ---- this is the generic plot file type name
        HeaderFile << versionName << '\n';

        HeaderFile << varnames.size() << '\n';

        for (int ivar = 0; ivar < varnames.size(); ++ivar) {
	    HeaderFile << varnames[ivar] << "\n";
        }
        HeaderFile << AMREX_SPACEDIM << '\n';
        HeaderFile << time << '\n';
        HeaderFile << finest_level << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << Geometry::ProbLo(i) << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << Geometry::ProbHi(i) << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i < finest_level; ++i) {
            HeaderFile << ref_ratio[i][0] << ' ';
	}
        HeaderFile << '\n';
	for (int i = 0; i <= finest_level; ++i) {
	    HeaderFile << geom[i].Domain() << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << level_steps[i] << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            for (int k = 0; k < AMREX_SPACEDIM; ++k) {
                HeaderFile << geom[i].CellSize()[k] << ' ';
	    }
            HeaderFile << '\n';
        }
        HeaderFile << (int) Geometry::Coord() << '\n';
        HeaderFile << "0\n";

	for (int level = 0; level <= finest_level; ++level) {
	    HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
	    HeaderFile << level_steps[level] << '\n';

	    for (int i = 0; i < bArray[level].size(); ++i)
	    {
		const Box &b(bArray[level][i]);
		RealBox loc = RealBox(b, geom[level].CellSize(), geom[level].ProbLo());
		for (int n = 0; n < AMREX_SPACEDIM; ++n) {
		    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
		}
	    }

	    HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
	}
}


void
WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
                         const Vector<const MultiFab*>& mf,
                         const Vector<std::string>& varnames,
                         const Vector<Geometry>& geom, Real time, const Vector<int>& level_steps,
                         const Vector<IntVect>& ref_ratio,
                         const std::string &versionName,
                         const std::string &levelPrefix,
                         const std::string &mfPrefix,
                         const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    int finest_level = nlevels-1;

//    int saveNFiles(VisMF::GetNOutFiles());
//    VisMF::SetNOutFiles(std::max(1024,saveNFiles));

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
      VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
      std::string HeaderFileName(plotfilename + "/Header");
      std::ofstream HeaderFile;
      HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
      HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
	                                      std::ofstream::trunc |
                                              std::ofstream::binary);
      if( ! HeaderFile.good()) {
        FileOpenFailed(HeaderFileName);
      }

      Vector<BoxArray> boxArrays(nlevels);
      for(int level(0); level < boxArrays.size(); ++level) {
	boxArrays[level] = mf[level]->boxArray();
      }

      WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, varnames,
                                 geom, time, level_steps, ref_ratio, versionName, levelPrefix, mfPrefix);
    }


    for (int level = 0; level <= finest_level; ++level)
    {
        const MultiFab* data;
        std::unique_ptr<MultiFab> mf_tmp;
        if (mf[level]->nGrow() > 0) {
            mf_tmp.reset(new MultiFab(mf[level]->boxArray(),
                                      mf[level]->DistributionMap(),
                                      mf[level]->nComp(), 0, MFInfo(),
                                      mf[level]->Factory()));
            MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
            data = mf_tmp.get();
        } else {
            data = mf[level];
        }

	    VisMF::Write(*data, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

//    VisMF::SetNOutFiles(saveNFiles);
}

#if defined(BL_USE_MPI) && defined(BL_USE_PNETCDF)
void WriteMultiLevelPlotfile_PNC(const std::string &plotfilename, int nlevels,
                                 const Vector<const MultiFab *> &mf,
                                 const Vector<std::string> &varnames,
                                 const Vector<Geometry> &geom, Real time, const Vector<int> &level_steps,
                                 const Vector<IntVect> &ref_ratio,
                                 const std::string &versionName,
                                 const std::string &levelPrefix,
                                 const std::string &mfPrefix,
                                 const Vector<std::string> &extra_dirs) {
    int err;
    int rank;
    int i, j, k;
    int ncid;
    int dimids[2];
    int varids[5];
    int finest_level(nlevels - 1);
    Real plo[AMREX_SPACEDIM], phi[AMREX_SPACEDIM];
    MPI_Offset start[2], count[2];
    std::vector<int> rratio(finest_level);
    std::vector<int> lvl_size(nlevels);
    std::vector<int> lvl_offs(nlevels);
    std::vector<Real> csize(nlevels);
    char name[1024];
    std::string fname;

    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size() + 1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    rank = ParallelDescriptor::MyProc();
    fname = std::string(plotfilename);
    if ( !fname.empty() && fname.back() == '/') {
	    fname = fname.substr(0, fname.size() - 1);
    }

    err = ncmpi_create(MPI_COMM_WORLD, fname.c_str(), 0, MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_create fail");
    }

    Vector<BoxArray> boxArrays(nlevels);
    for (int level(0); level < boxArrays.size(); ++level)
    {
        boxArrays[level] = mf[level]->boxArray();
    }
    lvl_offs[0] = 0;
    for (i = 0; i <= finest_level; i++) {
        lvl_size[i] = boxArrays[i].size();
        lvl_offs[i + 1] = lvl_offs[i] + lvl_size[i];
	}
    std::vector<std::vector<Real>> lvl_box_lo(lvl_offs[nlevels]);
    std::vector<std::vector<Real>> lvl_box_hi(lvl_offs[nlevels]);
    for (i = 0; i <= finest_level; i++) {
	    for (j = 0; j < lvl_size[i]; j++) {
            const Box &b(boxArrays[i][j]);
		    RealBox loc = RealBox(b, geom[i].CellSize(), geom[i].ProbLo());
            
            lvl_box_lo[lvl_offs[i] + j] = std::vector<Real>(AMREX_SPACEDIM);
            lvl_box_hi[lvl_offs[i] + j] = std::vector<Real>(AMREX_SPACEDIM);
            for (int k = 0; k < AMREX_SPACEDIM; k++) {
                lvl_box_lo[lvl_offs[i] + j][k] = loc.lo(k);
                lvl_box_hi[lvl_offs[i] + j][k] = loc.hi(k);
            }
	    }
	}

    // Define dimensions
    err = ncmpi_def_dim(ncid, "nlevels", nlevels, dimids);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_def_dim fail");
    }

    err = ncmpi_def_dim(ncid, "AMREX_SPACEDIM", AMREX_SPACEDIM, dimids + 1);
    if (err == NC_ENAMEINUSE){
        err = ncmpi_inq_dimid(ncid, "AMREX_SPACEDIM", dimids + 1);
        if (err != NC_NOERR){
            amrex::Error("ncmpi_inq_dimid fail");
        }
    }
    else if (err != NC_NOERR){
        amrex::Error("ncmpi_def_dim fail");
    }
    
    err = ncmpi_def_dim(ncid, "lvl_size_total", lvl_offs[nlevels], dimids + 2);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_def_dim fail");
    }

    // Variables
    err = ncmpi_def_var(ncid, "geom_lo", NC_INT, 2, dimids, varids);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_def_var fail");
    }

    err = ncmpi_def_var(ncid, "geom_hi", NC_INT, 2, dimids, varids + 1);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_def_var fail");
    }

#ifdef BL_USE_FLOAT
    err = ncmpi_def_var(ncid, "geom_CellSize", NC_FLOAT, 2, dimids, varids + 2);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_def_var fail");
    }

    dimids[0] = dimids[2];

    err = ncmpi_def_var(ncid, "lvl_box_lo", NC_FLOAT, 2, dimids, varids + 3);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_def_var fail");
    }

    err = ncmpi_def_var(ncid, "lvl_box_hi", NC_FLOAT, 2, dimids, varids + 4);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_def_var fail");
    }
#else
    err = ncmpi_def_var(ncid, "geom_CellSize", NC_DOUBLE, 2, dimids, varids + 2);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_def_var fail");
    }

    dimids[0] = dimids[2];

    err = ncmpi_def_var(ncid, "lvl_box_lo", NC_DOUBLE, 2, dimids, varids + 3);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_def_var fail");
    }

    err = ncmpi_def_var(ncid, "lvl_box_hi", NC_DOUBLE, 2, dimids, varids + 4);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_def_var fail");
    }
#endif

    // Attributes
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "geom_lo_varid", NC_INT, 1, varids);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "geom_hi_varid", NC_INT, 1, varids + 1);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "geom_CellSize_varid", NC_INT, 1, varids + 2);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "lvl_box_lo_varid", NC_INT, 1, varids + 3);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "lvl_box_hi_varid", NC_INT, 1, varids + 4);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "versionName", versionName.size(), versionName.c_str());
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_text fail");
    }

    i = varnames.size();
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "nvar", NC_INT, 1, &i);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    for (i = 0; i < varnames.size(); ++i) {
        sprintf(name, "varname_%d", i);
        err = ncmpi_put_att_text(ncid, NC_GLOBAL, name, varnames[i].size(), varnames[i].c_str());
        if (err != NC_NOERR){
            amrex::Error("ncmpi_put_att_text fail");
        }
    }

    i = AMREX_SPACEDIM;
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "AMREX_SPACEDIM", NC_INT, 1, &i);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "finest_level", NC_INT, 1, &finest_level);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    i = Geometry::Coord();
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "Geometry_Coord", NC_INT, 1, &i);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    for (i = 0; i < AMREX_SPACEDIM; ++i) {
        plo[i] = Geometry::ProbLo(i);
        phi[i] = Geometry::ProbHi(i);
	}
#ifdef BL_USE_FLOAT
    err = ncmpi_put_att_float(ncid, NC_GLOBAL, "Geometry_ProbLo", NC_FLOAT, AMREX_SPACEDIM, plo);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_float fail");
    }
    err = ncmpi_put_att_float(ncid, NC_GLOBAL, "Geometry_ProbHi", NC_FLOAT, AMREX_SPACEDIM, phi);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_float fail");
    }

    err = ncmpi_put_att_float(ncid, NC_GLOBAL, "time", NC_FLOAT, 1, &time);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_float fail");
    }
#else
    err = ncmpi_put_att_double(ncid, NC_GLOBAL, "Geometry_ProbLo", NC_DOUBLE, AMREX_SPACEDIM, plo);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_double fail");
    }
    err = ncmpi_put_att_double(ncid, NC_GLOBAL, "Geometry_ProbHi", NC_DOUBLE, AMREX_SPACEDIM, phi);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_double fail");
    }

    err = ncmpi_put_att_double(ncid, NC_GLOBAL, "time", NC_DOUBLE, 1, &time);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_double fail");
    }
#endif

    for (i = 0; i < finest_level; ++i) {
        rratio[i] = ref_ratio[i][0];
    }
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "ref_ratio", NC_INT, finest_level, rratio.data());
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "level_steps", NC_INT, nlevels, level_steps.data());
    if (err != NC_NOERR){
        amrex::Error("ncmpi_put_att_int fail");
    }

    // Data
    if (rank == 0){       
        start[1] = 0;
        count[0] = 1;
        count[1] = AMREX_SPACEDIM;
        for (i = 0; i <= finest_level; ++i) {
            start[0] = i;

            err = ncmpi_iput_vara_int(ncid, varids[0], start, count, ((const Box)(geom[i].Domain())).loVect(), NULL);
            if (err != NC_NOERR){
                amrex::Error("ncmpi_iput_vara_int fail");
            }
            
            err = ncmpi_iput_vara_int(ncid, varids[1], start, count, ((const Box)(geom[i].Domain())).hiVect(), NULL);
            if (err != NC_NOERR){
                amrex::Error("ncmpi_iput_vara_int fail");
            }

#ifdef BL_USE_FLOAT
            err = ncmpi_iput_vara_float(ncid, varids[2], start, count, geom[i].CellSize(), NULL);
            if (err != NC_NOERR){
                amrex::Error("ncmpi_iput_vara_float fail");
            }
            
            for (j = 0; j < lvl_size[i]; j++) {
                start[0] = lvl_offs[i] + j;
                err = ncmpi_iput_vara_float(ncid, varids[3], start, count, lvl_box_lo[lvl_offs[i] + j].data(), NULL);
                if (err != NC_NOERR){
                    amrex::Error("ncmpi_iput_vara_float fail");
                }
                err = ncmpi_iput_vara_float(ncid, varids[4], start, count, lvl_box_hi[lvl_offs[i] + j].data(), NULL);
                if (err != NC_NOERR){
                    amrex::Error("ncmpi_iput_vara_float fail");
                }
            }
#else
            err = ncmpi_iput_vara_double(ncid, varids[2], start, count, geom[i].CellSize(), NULL);
            if (err != NC_NOERR){
                amrex::Error("ncmpi_iput_vara_double fail");
            }

            for (j = 0; j < lvl_size[i]; j++) {
                start[0] = lvl_offs[i] + j;
                err = ncmpi_iput_vara_double(ncid, varids[3], start, count, lvl_box_lo[lvl_offs[i] + j].data(), NULL);
                if (err != NC_NOERR){
                    amrex::Error("ncmpi_iput_vara_double fail");
                }
                err = ncmpi_iput_vara_double(ncid, varids[4], start, count, lvl_box_hi[lvl_offs[i] + j].data(), NULL);
                if (err != NC_NOERR){
                    amrex::Error("ncmpi_iput_vara_double fail");
                }
            }
#endif
        }
    }

    err = ncmpi_enddef(ncid);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_enddef fail");
    }

    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_wait_all fail");
    }

    // Fab data
    for (int level = 0; level <= finest_level; ++level)
    {
        const MultiFab *data;
        std::unique_ptr<MultiFab> mf_tmp;
        if (mf[level]->nGrow() > 0)
        {
            mf_tmp.reset(new MultiFab(mf[level]->boxArray(),
                                      mf[level]->DistributionMap(),
                                      mf[level]->nComp(), 0, MFInfo(),
                                      mf[level]->Factory()));
            MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
            data = mf_tmp.get();
        }
        else
        {
            data = mf[level];
        }

        sprintf(name, "%s_%s_%d_%s", fname.c_str(), levelPrefix.c_str(), level, mfPrefix.c_str());
        VisMF::Write_PNC(ncid, *data, std::string(name), false);
    }

    err = ncmpi_close(ncid);
    if (err != NC_NOERR){
        amrex::Error("ncmpi_close fail");
    }
}
#endif

// write a plotfile to disk given:
// -plotfile name
// -vector of MultiFabs
// -vector of Geometrys
// variable names are written as "Var0", "Var1", etc.    
// refinement ratio is computed from the Geometry vector
// "time" and "level_steps" are set to zero
void WriteMLMF (const std::string &plotfilename,
                const Vector<const MultiFab*>& mf,
                const Vector<Geometry> &geom)
{
    int nlevs = mf.size();
    int ncomp = mf[0]->nComp();

    // variables names are "Var0", "Var1", etc.
    Vector<std::string> varnames(ncomp);
    for (int i=0; i<ncomp; ++i) {
        varnames[i] = "Var" + std::to_string(i);
    }

    // compute refinement ratio by looking at hi coordinates of domain at each level from
    // the geometry object
    Vector<IntVect> ref_ratio(nlevs-1);
    for (int i = 0; i < nlevs-1; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            int rr = (geom[i+1].Domain()).bigEnd(d)/(geom[i].Domain()).bigEnd(d);
            ref_ratio[i][d] = rr;
        }
    }

    // set step_array to zero
    Vector<int> step_array(nlevs,0);

    // set time to zero
    Real time = 0.;
    
    WriteMultiLevelPlotfile(plotfilename, nlevs, mf, varnames,
                            geom, time, step_array, ref_ratio);   
    
}    


void
WriteMultiLevelPlotfileHeaders (const std::string & plotfilename, int nlevels,
                                const Vector<const MultiFab*> & mf,
                                const Vector<std::string>     & varnames,
                                const Vector<Geometry>        & geom,
                                Real time, const Vector<int>  & level_steps,
                                const Vector<IntVect>     & ref_ratio,
                                const std::string         & versionName,
                                const std::string         & levelPrefix,
                                const std::string         & mfPrefix,
                                const Vector<std::string> & extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    int finest_level = nlevels-1;

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::string HeaderFileName(plotfilename + "/Header");
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                        std::ofstream::trunc |
                        std::ofstream::binary);
        if( ! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, varnames,
                                   geom, time, level_steps, ref_ratio, versionName, levelPrefix, mfPrefix);
    }

    for (int level = 0; level <= finest_level; ++level) {
        const MultiFab * data;
        data = mf[level];
        VisMF::WriteOnlyHeader(*data, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

}



void
WriteSingleLevelPlotfile (const std::string& plotfilename,
                          const MultiFab& mf, const Vector<std::string>& varnames,
                          const Geometry& geom, Real time, int level_step,
                          const std::string &versionName,
                          const std::string &levelPrefix,
                          const std::string &mfPrefix,
                          const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    WriteMultiLevelPlotfile(plotfilename, 1, mfarr, varnames, geomarr, time,
                            level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}

#ifdef AMREX_USE_EB
void
EB_WriteSingleLevelPlotfile (const std::string& plotfilename,
                             const MultiFab& mf, const Vector<std::string>& varnames,
                             const Geometry& geom, Real time, int level_step,
                             const std::string &versionName,
                             const std::string &levelPrefix,
                             const std::string &mfPrefix,
                             const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    EB_WriteMultiLevelPlotfile(plotfilename, 1, mfarr, varnames, geomarr, time,
                               level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}

void
EB_WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
                            const Vector<const MultiFab*>& mf,
                            const Vector<std::string>& varnames,
                            const Vector<Geometry>& geom, Real time, const Vector<int>& level_steps,
                            const Vector<IntVect>& ref_ratio,
                            const std::string &versionName,
                            const std::string &levelPrefix,
                            const std::string &mfPrefix,
                            const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf[0]->hasEBFabFactory(),
                                     "EB_WriteMultiLevelPlotfile: does not have EB Factory");

    int finest_level = nlevels-1;

//    int saveNFiles(VisMF::GetNOutFiles());
//    VisMF::SetNOutFiles(std::max(1024,saveNFiles));

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::string HeaderFileName(plotfilename + "/Header");
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
	                                        std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        Vector<std::string> vn = varnames;
        vn.push_back("vfrac");
        WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, vn,
                                   geom, time, level_steps, ref_ratio, versionName,
                                   levelPrefix, mfPrefix);

        for (int lev = 0; lev < nlevels; ++lev) {
            HeaderFile << "1.0e-6\n";
        }
    }


    for (int level = 0; level <= finest_level; ++level)
    {
        const int nc = mf[level]->nComp();
        MultiFab mf_tmp(mf[level]->boxArray(),
                        mf[level]->DistributionMap(),
                        nc+1, 0);
        MultiFab::Copy(mf_tmp, *mf[level], 0, 0, nc, 0);
        auto const& factory = dynamic_cast<EBFArrayBoxFactory const&>(mf[level]->Factory());
        MultiFab::Copy(mf_tmp, factory.getVolFrac(), 0, nc, 1, 0);
	VisMF::Write(mf_tmp, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

//    VisMF::SetNOutFiles(saveNFiles);
}

#endif
}

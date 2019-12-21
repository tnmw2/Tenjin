#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void PrintAllVarsTo1DGnuplotFile(CellArray &U, int picture, std::__cxx11::string filename, const Real *dx, const Real *prob_lo);


#ifdef AMREX_PARTICLES
std::unique_ptr<AmrTracerParticleContainer> AmrLevelAdv::TracerPC =  nullptr;
int AmrLevelAdv::do_tracers                       =  0;
#endif

//
//Default constructor.  Builds invalid object.
//
AmrLevelAdv::AmrLevelAdv ()
{
    flux_reg = 0;
}

//
//The basic constructor.
//
AmrLevelAdv::AmrLevelAdv (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time) 
{
    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,parameters.Ncomp);
        
}

//
//The destructor.
//
AmrLevelAdv::~AmrLevelAdv () 
{
    delete flux_reg;
}

//
//Restart from a checkpoint file.
//
void
AmrLevelAdv::restart (Amr&          papa,
	              std::istream& is,
                      bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,parameters.Ncomp);
}

void 
AmrLevelAdv::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old) 
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
#ifdef AMREX_PARTICLES
  if (do_tracers and level == 0) {
    TracerPC->Checkpoint(dir, "Tracer", true);
  }
#endif
}

//
//Write a plotfile to specified directory.
//
void AmrLevelAdv::writePlotFile (const std::string& dir, std::ostream& os, VisMF::How how)
{
    static int picture = 0;

    AmrLevel::writePlotFile (dir,os,how);

    /*MultiFab& S_new = get_new_data(Phi_Type);

    CellArray U(S_new,accessPattern,parameters);

    std::string filename = "Gnuplot/Barton1_REF";

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    PrintAllVarsTo1DGnuplotFile(U,picture,filename,dx,prob_lo);

    picture++;*/

#ifdef AMREX_PARTICLES
    if (do_tracers and level == 0) {
      TracerPC->Checkpoint(dir, "Tracer", true);
    }
#endif
}


//
//Cleanup data descriptors at end of run.
//
void
AmrLevelAdv::variableCleanUp () 
{
    desc_lst.clear();
#ifdef AMREX_PARTICLES
    TracerPC.reset();
#endif
}


//
//Initialize data on this level from another AmrLevelAdv (during regrid).
//
void
AmrLevelAdv::init (AmrLevel &old)
{
    AmrLevelAdv* oldlev = (AmrLevelAdv*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[Phi_Type].curTime();
    Real prev_time = oldlev->state[Phi_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(Phi_Type);

    FillPatch(old, S_new, 0, cur_time, Phi_Type, 0, parameters.Ncomp);
}

//
//Initialize data on this level after regridding if old level did not previously exist
//
void
AmrLevelAdv::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, parameters.Ncomp);
}

//
//Compute initial time step.
//
Real
AmrLevelAdv::initialTimeStep ()
{
    return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
AmrLevelAdv::computeInitialDt (int                   finest_level,
	  	               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Compute new `dt'.
//
void
AmrLevelAdv::computeNewDt (int                   finest_level,
		           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {
        AmrLevelAdv& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1) 
    {
	//
	// Limit dt's by pre-regrid dt
	//
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],dt_level[i]);
	}
    }
    else 
    {
	//
	// Limit dt's by change_max * old dt
	//
	static Real change_max = 1.1;
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
	}
    }
    
    //
    // Find the minimum over all levels
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Do work after timestep().
//
void
AmrLevelAdv::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();

#ifdef AMREX_PARTICLES    
    if (TracerPC)
      {
        const int ncycle = parent->nCycle(level);
	
        if (iteration < ncycle || level == 0)
	  {
            int ngrow = (level == 0) ? 0 : iteration;
	    
	    TracerPC->Redistribute(level, TracerPC->finestLevel(), ngrow);
	  }
      }
#endif
}

//
//Do work after regrid().
//
void
AmrLevelAdv::post_regrid (int lbase, int new_finest) {
#ifdef AMREX_PARTICLES
  if (TracerPC && level == lbase) {
      TracerPC->Redistribute(lbase);
  }
#endif
}

//
//Do work after a restart().
//
void
AmrLevelAdv::post_restart() 
{
#ifdef AMREX_PARTICLES
    if (do_tracers and level == 0) {
      BL_ASSERT(TracerPC == 0);
      TracerPC.reset(new AmrTracerParticleContainer(parent));
      TracerPC->Restart(parent->theRestartFile(), "Tracer");
    }
#endif
}

//
//Do work after init().
//
void
AmrLevelAdv::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

void
AmrLevelAdv::reflux ()
{
    BL_ASSERT(level<parent->finestLevel());

    const Real strt = amrex::second();

    getFluxReg(level+1).Reflux(get_new_data(Phi_Type),1.0,0,0,parameters.Ncomp,geom);
    
    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = amrex::second() - strt;
	
        ParallelDescriptor::ReduceRealMax(end,IOProc);
	
        amrex::Print() << "AmrLevelAdv::reflux() at level " << level 
                       << " : time = " << end << std::endl;
    }
}

void
AmrLevelAdv::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(Phi_Type);
}

void
AmrLevelAdv::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    AmrLevelAdv& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);
    
    amrex::average_down(S_fine,S_crse,
                         fine_lev.geom,geom,
                         0,S_fine.nComp(),parent->refRatio(level));
}

#ifdef AMREX_PARTICLES
void
AmrLevelAdv::init_particles ()
{
  if (do_tracers and level == 0)
    {
      BL_ASSERT(TracerPC == nullptr);
      
      TracerPC.reset(new AmrTracerParticleContainer(parent));
      TracerPC->do_tiling = true;
      TracerPC->tile_size = IntVect(AMREX_D_DECL(1024000,4,4));

      AmrTracerParticleContainer::ParticleInitData pdata = {AMREX_D_DECL(0.0, 0.0, 0.0)};

      TracerPC->SetVerbose(0);
      TracerPC->InitOnePerCell(0.5, 0.5, 0.5, pdata);

      TracerPC->Redistribute();
    }
}
#endif

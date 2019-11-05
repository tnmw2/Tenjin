#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

#include "simulationheader.h"


using namespace amrex;

int      AmrLevelAdv::verbose           = 0;
Real     AmrLevelAdv::cfl               = 0.9;
int      AmrLevelAdv::do_reflux         = 1;
Real	 AmrLevelAdv::advectionVeloctiy = 0.0;

ParameterStruct AmrLevelAdv::parameters;
InitialStruct   AmrLevelAdv::initial;
PlasticEOS      AmrLevelAdv::plastic;
AccessPattern 	AmrLevelAdv::accessPattern(AmrLevelAdv::parameters);

int      AmrLevelAdv::NUM_STATE       = 1;  // One variable in the state
int      AmrLevelAdv::NUM_GROW        = 1;  // number of ghost cells

void calc_5Wave_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& ULStarbox, BoxAccessCellArray& URStarbox, BoxAccessCellArray& UStarStarbox, ParameterStruct& parameters, Direction_enum d);
void calc_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters, Direction_enum d);
void update(BoxAccessCellArray& fluxbox, BoxAccessCellArray& Ubox, BoxAccessCellArray& U1box, ParameterStruct& parameters, Direction_enum d, Real dt, const Real* dx);
void MUSCLextrapolate(BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, BoxAccessCellArray& grad, Direction_enum d);


void C_nullfill(double* adv, const int& adv_lo, const int& adv_hi, const int& domlo, const int& domhi, const int* delta, const int* xlo, const double* time, const double*, const double*, const int* bc)
{

	return;
}

void AmrLevelAdv::initData ()
{

    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(Phi_Type);
    Real cur_time   = state[Phi_Type].curTime();

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }
     
    CellArray U(S_new,accessPattern,parameters); 
     
    setInitialConditions(U,parameters,initial,dx,prob_lo);
    

#ifdef AMREX_PARTICLES
    init_particles();
#endif

    if (verbose) {
	amrex::Print() << "Done initializing the level " << level 
                       << " data " << std::endl;
    }
}

void ScaleFluxes(Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Real dt, const Real* dx, AccessPattern& accessPattern, ParameterStruct& parameters)
{
	CellArray F0(flux_arr[0],accessPattern,parameters);
	CellArray F1(flux_arr[1],accessPattern,parameters);
	CellArray F2(flux_arr[2],accessPattern,parameters);
	
	for(MFIter mfi(flux_arr[0]); mfi.isValid(); ++mfi )
	{
		const Box& bx = mfi.validbox();

		BoxAccessCellArray  xFlux(mfi,bx,F0);
		
		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		for(auto n : accessPattern.conservativeVariables)
		{
			for 		(int k = lo.z; k <= hi.z; ++k)
			{
				for 	(int j = lo.y; j <= hi.y; ++j)
				{
					for (int i = lo.x; i <= hi.x; ++i)
					{
						xFlux(i,j,k,n) *= (dt*dx[1]);	
					}
				}
			}
		}
	}
	
	for(MFIter mfi(flux_arr[1]); mfi.isValid(); ++mfi )
	{
		const Box& bx = mfi.validbox();

		BoxAccessCellArray  yFlux(mfi,bx,F1);
		
		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		for(auto n : accessPattern.conservativeVariables)
		{
			for 		(int k = lo.z; k <= hi.z; ++k)
			{
				for 	(int j = lo.y; j <= hi.y; ++j)
				{
					for (int i = lo.x; i <= hi.x; ++i)
					{
						yFlux(i,j,k,n) *= (dt*dx[0]);	
					}
				}
			}
		}
	}
		
		
}

void AmrLevelAdv::AMR_HLLCadvance(MultiFab& S_new,CellArray& U,CellArray& U1, CellArray& UL, CellArray& UR, CellArray& MUSCLgrad, CellArray& ULStar, CellArray& URStar, CellArray& UStarStar, Array<MultiFab, AMREX_SPACEDIM>& flux_arr,ParameterStruct& parameters, const Real* dx, Real dt, Real time)
{
    Direction_enum d;
    
    U.conservativeToPrimitive();

    U1 = U;

    for(int dir = 0; dir < AMREX_SPACEDIM ; dir++)
    {
        d = (Direction_enum)dir;

        UL = U;
        UR = U;

        /*-------------------------------------------------------------
         * Perform MUSCL extrapolation.
         * -----------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(U.data); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            /*-------------------------------------------------------------
             * Data can't be accessed straight from a Multifab so we make
             * some wrappers to hold the FArrayBoxes that can access the
             * data called BoxAccessCellArray.
             * -----------------------------------------------------------*/

            BoxAccessCellArray  Ubox(mfi,bx,U);
            BoxAccessCellArray  ULbox(mfi,bx,UL);
            BoxAccessCellArray  URbox(mfi,bx,UR);


            if(parameters.MUSCL)
            {
                BoxAccessCellArray  gradbox(mfi,bx,MUSCLgrad);

                MUSCLextrapolate(Ubox,ULbox,URbox,gradbox,d);
            }

            /*if(parameters.THINC)
            {
                BoxAccessCellArray ULTHINC(mfi,bx,ULStar);
                BoxAccessCellArray URTHINC(mfi,bx,URStar);
                BoxAccessTHINCArray THINCbox(mfi,bx,THINC);

                THINCbox.THINCreconstruction(Ubox,ULbox,URbox,ULTHINC,URTHINC,parameters,d);

                UL.cleanUpAlpha();
                UR.cleanUpAlpha();

            }*/


            ULbox.primitiveToConservative();
            URbox.primitiveToConservative();

            ULbox.getSoundSpeed();
            URbox.getSoundSpeed();

        }

        /*FillDomainBoundary(UL.data, geom, bc);
        FillDomainBoundary(UR.data, geom, bc);

        UL.data.FillBoundary(geom.periodicity());
        UR.data.FillBoundary(geom.periodicity());*/
        
        FillPatch(*this, UL.data, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
		FillPatch(*this, UR.data, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
		


        /*-------------------------------------------------------------
         * Calulate HLLC flux and update the new array.
         * -----------------------------------------------------------*/

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(S_new); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            FArrayBox& flux_fab   = flux_arr[d][mfi];

            BoxAccessCellArray Ubox(mfi,bx,U);
            BoxAccessCellArray U1box(mfi,bx,U1);
            BoxAccessCellArray ULbox(mfi,bx,UL);
            BoxAccessCellArray URbox(mfi,bx,UR);
            BoxAccessCellArray ULStarbox(mfi,bx,ULStar);
            BoxAccessCellArray URStarbox(mfi,bx,URStar);
            BoxAccessCellArray UStarStarbox(mfi,bx,UStarStar);
            BoxAccessCellArray fluxbox(bx,flux_fab,U); 

            if(parameters.SOLID)
            {
                calc_5Wave_fluxes(fluxbox, ULbox, URbox, ULStarbox, URStarbox, UStarStarbox, parameters,d);
            }
            else
            {
                calc_fluxes(fluxbox, ULbox, URbox, ULStarbox, parameters,d);
            }

            update(fluxbox, Ubox, U1box, parameters,d,dt,dx);

        }
    }

    U1.conservativeToPrimitive();
    /*FillDomainBoundary(U1.data, geom, bc);
    U1.data.FillBoundary(geom.periodicity());*/

    //FillPatch(*this, U1.data, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);


}

Real AmrLevelAdv::advance (Real time, Real dt, int  iteration, int  ncycle)
{

    MultiFab& S_new = get_new_data(Phi_Type);

    const Real prev_time = state[Phi_Type].prevTime();
    const Real cur_time = state[Phi_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    //Print() << "In advance: " << level << " " << dx[0] << std::endl;


    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;
    
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
    {
		fine = &getFluxReg(level+1);
		fine->setVal(0.0);
    }

    if (do_reflux && level > 0)
    {
		current = &getFluxReg(level);
    }

    Array <MultiFab, AMREX_SPACEDIM> fluxes;
    Array <MultiFab, AMREX_SPACEDIM> fluxes1;
    Array <MultiFab, AMREX_SPACEDIM> fluxes2;



    //if(do_reflux)
    //{
		for (int j = 0; j < BL_SPACEDIM; j++)
		{
			BoxArray ba = S_new.boxArray();
			ba.surroundingNodes(j);
            fluxes [j].define(ba, dmap, parameters.Ncomp, 0);
            fluxes [j].setVal(0.0);
            fluxes1[j].define(ba, dmap, parameters.Ncomp, 0);
            fluxes1[j].setVal(0.0);
            fluxes2[j].define(ba, dmap, parameters.Ncomp, 0);
            fluxes2[j].setVal(0.0);
		}
    //}

    // States with ghost cells
    MultiFab S0(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, S0, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    MultiFab S1(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, S1, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    MultiFab S2(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, S2, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    MultiFab SL(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, SL, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    MultiFab SR(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, SR, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    MultiFab SLStar(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, SLStar, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    MultiFab SRStar(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, SRStar, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    MultiFab SStarStar(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, SStarStar, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    MultiFab Sgrad(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, Sgrad, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
         
    
    CellArray U (S0,accessPattern,parameters);
    CellArray U1(S1,accessPattern,parameters);
    CellArray U2(S2,accessPattern,parameters);
    CellArray UL(SL,accessPattern,parameters);
    CellArray UR(SR,accessPattern,parameters);
    CellArray ULStar(SLStar,accessPattern,parameters);
    CellArray URStar(SRStar,accessPattern,parameters);
	CellArray UStarStar(SStarStar,accessPattern,parameters);
    CellArray MUSCLgrad(Sgrad,accessPattern,parameters);
    
    AMR_HLLCadvance(S_new,U,U1,UL,UR,MUSCLgrad,ULStar,URStar,UStarStar,fluxes1,parameters,dx,dt,time);
    AMR_HLLCadvance(S_new,U1,U2,UL,UR,MUSCLgrad,ULStar,URStar,UStarStar,fluxes2,parameters,dx,dt,time);

    U1 = ((U*(1.0/2.0))+(U2*(1.0/2.0)));

    for (int j = 0; j < BL_SPACEDIM; j++)
    {
        MultiFab::LinComb(fluxes[j],0.5,fluxes1[j],0,0.5,fluxes2[j],0,0,fluxes[j].nComp(),0);
    }

    ScaleFluxes(fluxes,dt,dx,accessPattern,parameters);

    MultiFab::Copy(S_new, U1.data, 0, 0, U1.data.nComp(), S_new.nGrow());

/*#ifdef _OPENMP
#pragma omp parallel
#endif
    {
		//FArrayBox flux[BL_SPACEDIM];

			for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
			{
				const Box& bx = mfi.tilebox();
				
				
				

				/*const FArrayBox& statein  = Sborder[mfi];
				FArrayBox&       stateout =   S_new[mfi];
				
				Array4<Real const> const& Uold = statein.array();
				Array4<Real> const& Unew = stateout.array();
				
				// Allocate fabs for fluxes and Godunov velocities.
				for (int i = 0; i < BL_SPACEDIM ; i++)
				{
					const Box& bxtmp = amrex::surroundingNodes(bx,i);
					flux[i].resize(bxtmp,parameters.Ncomp);
				}*/
				
				/*Array4<Real> const& xflux = flux[0].array();
				Array4<Real> const& yflux = flux[1].array();
				Array4<Real> const& zflux = flux[2].array();
  
				C_advect(time,bx,Uold,Unew,xflux,yflux,zflux,dx,dt,cfl,parameters.Ncomp,advectionVeloctiy);	
			

			if (do_reflux)
			{
				for (int i = 0; i < BL_SPACEDIM ; i++)
				{
					fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));
				}
			}
		}
    }*/
    
    if (do_reflux)
    {
		if (current)
		{
			for (int i = 0; i < BL_SPACEDIM ; i++)
			{
				current->FineAdd(fluxes[i],i,0,0,parameters.Ncomp,1.);
			}
		}
		if (fine)
		{
			for (int i = 0; i < BL_SPACEDIM ; i++)
			{
				fine->CrseInit(fluxes[i],i,0,0,parameters.Ncomp,-1.);
			}
		}
    }

#ifdef AMREX_PARTICLES
    if (TracerPC) {
      TracerPC->AdvectWithUmac(Umac, level, dt);
    }
#endif

    return dt;
}

Real AmrLevelAdv::estTimeStep (Real)
{
    // This is just a dummy value to start with 
    Real dt_est  = 1.0e+20;

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real cur_time = state[Phi_Type].curTime();
    
    MultiFab& S_new = get_new_data(Phi_Type);
    
    CellArray U(S_new,accessPattern,parameters);

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
		for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();

            BoxAccessCellArray  Ubox(mfi,bx,U);
            
			const auto lo = lbound(bx);
			const auto hi = ubound(bx);
			
			Ubox.getSoundSpeed();			

			for 		(int k = lo.z; k <= hi.z; ++k)
			{
				for 	(int j = lo.y; j <= hi.y; ++j)
				{
					for (int i = lo.x; i <= hi.x; ++i)
					{
						for(int n=0;n<AMREX_SPACEDIM;n++)
						{
							dt_est = std::min(dt_est, dx[n] / (Ubox(i,j,k,SOUNDSPEED) + fabs(Ubox(i,j,k,VELOCITY,0,n)))); //(advectionVeloctiy));//
						}
					}
				}
			}
		}
    }

    ParallelDescriptor::ReduceRealMin(dt_est);
    dt_est *= cfl;

    if (verbose) {
	amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level 
                       << ":  dt_est = " << dt_est << std::endl;
    }
    
    return dt_est;
}

Real findBiggestGradientOnLevel(const Box& box, BoxAccessCellArray& U, const Real* dx, Real time, Real level)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Real ax,ay,gradient;
    
    Real max = 0.0;
    
    for 		(int k = lo.z; k <= hi.z; ++k)
	{
		for 	(int j = lo.y; j <= hi.y; ++j)
		{
			for (int i = lo.x; i <= hi.x; ++i)
			{
				
				ax = (U(i+1,j,k,RHO)-U(i-1,j,k,RHO))/(2.0*dx[0]);
				ay = (U(i,j+1,k,RHO)-U(i,j-1,k,RHO))/(2.0*dx[1]);

				gradient = sqrt(ax*ax+ay*ay);
				
				max = (gradient > max ? gradient : max);
			}
		}
	}
	
	return max;			
}

void C_state_error(Array4<char> const& tagarr, const Box& box, BoxAccessCellArray& U, const Real* dx, Real time, Real level, Real gradMax)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    

    Vector<Real> levelSpecificCoefficient{0.3,0.3,0.3,0.3,0.3,0.3};
    
    Real ax,ay,az,gradient;

    for 		(int k = lo.z; k <= hi.z; ++k)
	{
		for 	(int j = lo.y; j <= hi.y; ++j)
		{
			for (int i = lo.x; i <= hi.x; ++i)
			{
				
				ax = (U(i+1,j,k,RHO)-U(i-1,j,k,RHO))/(2.0*dx[0]);
				ay = (U(i,j+1,k,RHO)-U(i,j-1,k,RHO))/(2.0*dx[1]);

				gradient = sqrt(ax*ax+ay*ay);
				

                if(gradient >= levelSpecificCoefficient[level]*gradMax)
				{
					tagarr(i,j,k) = TagBox::SET;
                }
			}	
		}
	}

	return;

}

void AmrLevelAdv::errorEst (TagBoxArray& tags, int clearval, int tagval, Real time, int n_error_buf, int ngrow)
{
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(Phi_Type);
    
    // State with ghost cells
    MultiFab Sborder(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    
    CellArray U(Sborder,accessPattern,parameters);
    
    Real gradMax = 0.0;
    Real temp    = 0.0;


#ifdef _OPENMP
#pragma omp parallel reduction(max:gradMax)
#endif
    {
		for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
		{
			const Box& bx	  			= mfi.validbox();
			
			BoxAccessCellArray 			Ubox(mfi,bx,U);

            temp = findBiggestGradientOnLevel(bx,Ubox,dx,time,level);
            
            gradMax = ( (temp > gradMax) ? temp : gradMax );
		}
	}
	
	ParallelDescriptor::ReduceRealMax(gradMax);
	
    Print() << level << " " << gradMax << std::endl;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
		for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
		{
			const Box& bx	  			= mfi.validbox();
			
			BoxAccessCellArray 			Ubox(mfi,bx,U);

            TagBox&    tagfab	  		= tags[mfi];
            
            Array4<char> const& tagarr 	= tagfab.array();
            
            C_state_error(tagarr,bx,Ubox,dx,time,level,gradMax);
            
        }
	}
}

void FStyle_State_error(Array4<char> const& tag, const int* tag_lo, const int* tag_hi, BoxAccessCellArray& state, const Box& box, int set, int clear, const Real* dx, const Real* prob_lo, Real time, int level)
{
	int dim;
	
	Real ax,ay,az,gradient;
	
	const auto lo = lbound(box);
	const auto hi = ubound(box);
	
	if(lo.z == hi.z)
	{
		dim = 2;
	}
	else
	{
		dim = 3;
	}
	
    if(level < 10)
	{
		for 		(int k = lo.z; k <= hi.z; ++k)
		{
			for 	(int j = lo.y; j <= hi.y; ++j)
			{
				for (int i = lo.x; i <= hi.x; ++i)
				{
					
					ax = (state(i+1,j,k,RHO)-state(i-1,j,k,RHO))/(2.0*dx[0]);
					ay = (state(i,j+1,k,RHO)-state(i,j-1,k,RHO))/(2.0*dx[1]);

					
					if(dim == 2)
					{
						az = 0.0;
					}
					else
					{
						Print() << "No tagging in 3D yet" << std::endl;
						exit(1);
					}
					
					gradient = sqrt(ax*ax+ay*ay+az*az);
					

					if(gradient >= 2.0)
					{
						tag(i,j,k) = TagBox::SET;		
					}
					else
					{
						tag(i,j,k) = TagBox::CLEAR;
					}
				}
			}
		}
	}
}


/*void
AmrLevelAdv::errorEst (TagBoxArray& tags, int clearval, int tagval, Real time, int n_error_buf, int ngrow)
{
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(Phi_Type);
    
    // State with ghost cells
    MultiFab Sborder(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    CellArray U(Sborder,accessPattern,parameters);
    
    
    

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;
	
		for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
		{
			const Box&  tilebx  = mfi.tilebox();

			TagBox&     tagfab  = tags[mfi];
			
			BoxAccessCellArray Ubox(mfi,tilebx,U);
			
			Array4<char> const& tagarr 	= tagfab.array();

			int*        tptr    = itags.dataPtr();
			
			const int*  tlo     = tilebx.loVect();
			const int*  thi     = tilebx.hiVect();
			
			FStyle_State_error(tagarr,tlo,thi,Ubox,tilebx,tagval,clearval,dx,prob_lo,time,level);
		}
    }
}*/

void AmrLevelAdv::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    read_params();

    desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),StateDescriptor::Point,0,parameters.Ncomp,&cell_cons_interp);

    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    
    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
		lo_bc[i] = hi_bc[i] = BCType::foextrap;   // transmissive boundaries
    }
    
    BCRec bc(lo_bc, hi_bc);
    
    for(int n = 0; n<parameters.Ncomp; n++)
    {
		desc_lst.setComponent(Phi_Type, n, accessPattern.variableNames[n] , bc, 
			  StateDescriptor::BndryFunc(phifill)); 
	}
}

void AmrLevelAdv::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("adv");   

    pp.query("v",verbose);
    pp.query("cfl",cfl);
    pp.query("do_reflux",do_reflux);
    pp.query("vel",advectionVeloctiy);
    
    libConfigInitialiseDataStructs(parameters,initial,plastic);
    
    accessPattern.define(parameters);
    
    parameters.Ncomp = accessPattern.variableNames.size();

    Geometry const* gg = AMReX::top()->getDefaultGeometry();

    // This tutorial code only supports Cartesian coordinates.
    if (! gg->IsCartesian()) {
	amrex::Abort("Please set geom.coord_sys = 0");
    }
}

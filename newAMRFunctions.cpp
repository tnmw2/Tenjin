#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

#include "simulationheader.h"

using namespace amrex;

/*-------------------------------------------------------------
 * Declaring static member data
 * -----------------------------------------------------------*/

int      AmrLevelAdv::verbose           = 0;
Real     AmrLevelAdv::cfl               = 0.9;
int      AmrLevelAdv::do_reflux         = 1;
Real	 AmrLevelAdv::advectionVeloctiy = 0.0;

ParameterStruct AmrLevelAdv::parameters;
InitialStruct   AmrLevelAdv::initial;
PlasticEOS      AmrLevelAdv::plastic;
AccessPattern 	AmrLevelAdv::accessPattern(AmrLevelAdv::parameters);
Vector<Real>    AmrLevelAdv::levelGradientCoefficients;
Vector<BCRec>   AmrLevelAdv::bc;

int      AmrLevelAdv::NUM_STATE       = 1;  // One variable in the state
int      AmrLevelAdv::NUM_GROW        = 2;  // number of ghost cells

/*-------------------------------------------------------------
 * Non-AMR functions in other files that can interface to
 * the new AMR_HLLC functions
 * -----------------------------------------------------------*/

void calc_5Wave_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& ULStarbox, BoxAccessCellArray& URStarbox, BoxAccessCellArray& UStarStarbox, ParameterStruct& parameters, Direction_enum d,const Real* dx, const Real* prob_lo);
void calc_fluxes(BoxAccessCellArray& fluxbox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters, Direction_enum d,const Real* dx, const Real* prob_lo);
void update(BoxAccessCellArray& fluxbox, BoxAccessCellArray& Ubox, BoxAccessCellArray& U1box, ParameterStruct& parameters, Direction_enum d, Real dt, const Real* dx);
void MUSCLextrapolate(const BoxAccessCellArray& U, BoxAccessCellArray& UL, BoxAccessCellArray& UR, Direction_enum d);

void customAbort(Vector<Real>& values, std::string& Message)
{
    std::ostringstream stream;

    stream << std::endl;

    for(auto n : values)
    {
        stream << n << std::endl;
    }

    std::string error = Message + stream.str();

    Abort(error);
}

void smooth(BoxAccessCellArray& U, BoxAccessCellArray& U1, ParameterStruct& parameters, InitialStruct& initial,const Real* dx, const Real* prob_lo)
{
    const auto lo = lbound(U.box);
    const auto hi = ubound(U.box);

    Real sum     = 0.0;
    int  counter = 0;

    Vector<int> pm {-1,0,1};

    Vector<MaterialSpecifier> smoothing;

    smoothing.push_back(MaterialSpecifier (ALPHA,0,0,0));
    smoothing.push_back(MaterialSpecifier (ALPHA,1,0,0));
    smoothing.push_back(MaterialSpecifier (ALPHA,2,0,0));

    for(auto n : smoothing)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    Real sum     = 0.0;
                    int  counter = 0;


                    for(auto nj : pm)
                    {
                        for(auto ni : pm)
                        {
                            sum += U.neighbour(ni,nj,0,i,j,k,n);

                            counter++;
                        }
                    }

                    /*if(std::abs(((1.0/counter)*(sum) - U1(i,j,k,n))/U1(i,j,k,n)) > 0.1 )
                    {
                        Print() << (1.0/counter)*(sum) << " " << U1(i,j,k,n) << std::endl;
                    }*/

                    U1(i,j,k,n) = (1.0/counter)*(sum);
                }
            }
        }
    }
}

void AmrLevelAdv::initData ()
{

    const Real* dx      = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new     = get_new_data(Phi_Type);

    if (verbose)
    {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

    /*-------------------------------------------------------------
     * The old setInitialConditions function can still be used to
     * initialise data on a level, but now works in terms of real
     * position instead of array coordinates
     * -----------------------------------------------------------*/
     
    CellArray U(S_new,accessPattern,parameters); 


    setInitialConditions(U,parameters,initial,dx,prob_lo);

    #ifdef AMREX_PARTICLES
    init_particles();
    #endif

    if (verbose)
    {
        amrex::Print() << "Done initializing the level " << level  << " data " << std::endl;
    }
}

void ScaleFluxes(Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Real dt, const Real* dx, AccessPattern& accessPattern, ParameterStruct& parameters)
{
    /*-------------------------------------------------------------
     * Fluxes need to be rescaled for the refluxing operation.
     * We put CellArray wrappers around the flux mulitfabs
     * to access them easily
     * -----------------------------------------------------------*/

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
                        if(n.var == V_TENSOR || n.var == ALPHA)
                        {
                            xFlux(i,j,k,n) *= 0.0;
                        }
                        else
                        {
                            xFlux(i,j,k,n) *= (dt*dx[1]);
                        }
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
                        if(n.var == V_TENSOR || n.var == ALPHA)
                        {
                            yFlux(i,j,k,n) *= 0.0;
                        }
                        else
                        {
                            yFlux(i,j,k,n) *= (dt*dx[0]);
                        }
					}
				}
			}
		}
	}	
}

void AmrLevelAdv::SINGLE_HLLC_advance(MultiFab& S_new, Array<MultiFab, AMREX_SPACEDIM>& fluxes, ParameterStruct& parameters, const Real* dx, const Real* prob_lo, Real dt, Real time)
{
    Array <MultiFab, AMREX_SPACEDIM> fluxes1;
    Array <MultiFab, AMREX_SPACEDIM> fluxes2;

    for (int j = 0; j < BL_SPACEDIM; j++)
    {
        BoxArray ba = S_new.boxArray();
        ba.surroundingNodes(j);

        fluxes1[j].define(ba, dmap, parameters.Ncomp, 0);
        fluxes2[j].define(ba, dmap, parameters.Ncomp, 0);

    }


    const int TWOGHOST = 2;
    const int ONEGHOST = 1;

    MultiFab S0(grids, dmap, parameters.Ncomp, TWOGHOST);
    FillPatch(*this, S0, TWOGHOST, time, Phi_Type, 0, parameters.Ncomp);
    MultiFab S1(grids, dmap, parameters.Ncomp, TWOGHOST);

    MultiFab SL(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SR(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SLStar(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SRStar(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SStarStar(grids, dmap, parameters.Ncomp, ONEGHOST);

    THINCArray THINCArr(this->grids,this->dmap,ONEGHOST,parameters);

    CellArray U (S0,accessPattern,parameters);
    CellArray U1(S1,accessPattern,parameters);
    CellArray UL(SL,accessPattern,parameters);
    CellArray UR(SR,accessPattern,parameters);
    CellArray ULStar (SLStar,accessPattern,parameters);
    CellArray URStar (SRStar,accessPattern,parameters);
    CellArray UStarStar (SRStar,accessPattern,parameters);

    AMR_HLLCadvance(S_new,U,U1,UL,UR,ULStar,URStar,UStarStar,fluxes,fluxes1,fluxes2,THINCArr,parameters,dx,prob_lo,dt,time);

    MultiFab::Copy(S_new, U1.data, 0, 0, U1.data.nComp(), S_new.nGrow());
    FillPatch(*this, U1.data, TWOGHOST, time, Phi_Type, 0, U1.data.nComp());

}

void AmrLevelAdv::RK2_HLLC_advance(MultiFab& S_new, Array<MultiFab, AMREX_SPACEDIM>& fluxes, ParameterStruct& parameters, const Real* dx, const Real* prob_lo, Real dt, Real time)
{

    const int TWOGHOST = 2;
    const int ONEGHOST = 1;

    /*-------------------------------------------------------------
   * Declared flux multifabs. Several are needed as we need to
   * keep track of the flux in each Runge-Kutta timestep to such
   * that we can reflux appropriately
   * -----------------------------------------------------------*/

    Array <MultiFab, AMREX_SPACEDIM> fluxes1;
    Array <MultiFab, AMREX_SPACEDIM> fluxes2;

    for (int j = 0; j < BL_SPACEDIM; j++)
    {
        BoxArray ba = S_new.boxArray();
        ba.surroundingNodes(j);

        fluxes1[j].define(ba, dmap, parameters.Ncomp, 0);
        fluxes2[j].define(ba, dmap, parameters.Ncomp, 0);

    }


    MultiFab S0(grids, dmap, parameters.Ncomp, TWOGHOST);
    FillPatch(*this, S0, TWOGHOST, time, Phi_Type, 0, parameters.Ncomp);
    MultiFab S1(grids, dmap, parameters.Ncomp, TWOGHOST);
    MultiFab S2(grids, dmap, parameters.Ncomp, TWOGHOST);

    MultiFab SL(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SR(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SLStar(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SRStar(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SStarStar(grids, dmap, parameters.Ncomp, ONEGHOST);

    THINCArray THINCArr(this->grids,this->dmap,ONEGHOST,parameters);

    CellArray U (S0,accessPattern,parameters);
    CellArray U1(S1,accessPattern,parameters);
    CellArray U2(S2,accessPattern,parameters);
    CellArray UL(SL,accessPattern,parameters);
    CellArray UR(SR,accessPattern,parameters);
    CellArray ULStar (SLStar,accessPattern,parameters);
    CellArray URStar (SRStar,accessPattern,parameters);
    CellArray UStarStar (SRStar,accessPattern,parameters);

    AMR_HLLCadvance(S_new,U,U1,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes,fluxes2,THINCArr,parameters,dx,prob_lo,dt,time);

    AMR_HLLCadvance(S_new,U1,U2,UL,UR,ULStar,URStar,UStarStar,fluxes2, fluxes, fluxes1,THINCArr,parameters,dx,prob_lo,dt,time);

    U1 = ((U*(1.0/2.0))+(U2*(1.0/2.0)));

    MultiFab::Copy(S_new, U1.data, 0, 0, U1.data.nComp(), S_new.nGrow());
    FillPatch(*this, U1.data, TWOGHOST, time, Phi_Type, 0, U1.data.nComp());

    /*if(do_reflux)
    {
        for (int j = 0; j < BL_SPACEDIM; j++)
        {
            MultiFab::LinComb(fluxes[j],0.5,fluxes1[j],0,0.5,fluxes2[j],0,0,fluxes[j].nComp(),0);
        }

        ScaleFluxes(fluxes,dt,dx,accessPattern,parameters);
    }*/


}

void AmrLevelAdv::RK3_HLLC_advance(MultiFab& S_new, Array<MultiFab, AMREX_SPACEDIM>& fluxes, ParameterStruct& parameters, const Real* dx, const Real* prob_lo, Real dt, Real time)
{

    const int TWOGHOST = 2;
    const int ONEGHOST = 1;

    /*-------------------------------------------------------------
   * Declared flux multifabs. Several are needed as we need to
   * keep track of the flux in each Runge-Kutta timestep to such
   * that we can reflux appropriately
   * -----------------------------------------------------------*/

    Array <MultiFab, AMREX_SPACEDIM> fluxes1;
    Array <MultiFab, AMREX_SPACEDIM> fluxes2;
    Array <MultiFab, AMREX_SPACEDIM> fluxes3;

    for (int j = 0; j < BL_SPACEDIM; j++)
    {
        BoxArray ba = S_new.boxArray();
        ba.surroundingNodes(j);

        fluxes1[j].define(ba, dmap, parameters.Ncomp, 0);
        fluxes2[j].define(ba, dmap, parameters.Ncomp, 0);
        fluxes3[j].define(ba, dmap, parameters.Ncomp, 0);

    }


    MultiFab S0(grids, dmap, parameters.Ncomp, TWOGHOST);
    FillPatch(*this, S0, TWOGHOST, time, Phi_Type, 0, parameters.Ncomp);

    MultiFab S1(grids, dmap, parameters.Ncomp, TWOGHOST);
    MultiFab S2(grids, dmap, parameters.Ncomp, TWOGHOST);
    MultiFab S3(grids, dmap, parameters.Ncomp, TWOGHOST);

    MultiFab SL(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SR(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SLStar(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SRStar(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SStarStar(grids, dmap, parameters.Ncomp, ONEGHOST);

    THINCArray THINCArr(this->grids,this->dmap,ONEGHOST,parameters);

    CellArray U (S0,accessPattern,parameters);
    CellArray U1(S1,accessPattern,parameters);
    CellArray U2(S2,accessPattern,parameters);
    CellArray U3(S3,accessPattern,parameters);
    CellArray UL(SL,accessPattern,parameters);
    CellArray UR(SR,accessPattern,parameters);
    CellArray ULStar (SLStar,accessPattern,parameters);
    CellArray URStar (SRStar,accessPattern,parameters);
    CellArray UStarStar (SRStar,accessPattern,parameters);

    AMR_HLLCadvance(S_new,U,U1,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes2,fluxes3,THINCArr,parameters,dx,prob_lo,dt,time);

    AMR_HLLCadvance(S_new,U1,U2,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes2,fluxes3,THINCArr,parameters,dx,prob_lo,dt,time);

    U1 = ((U*(3.0/4.0))+(U2*(1.0/4.0)));

    AMR_HLLCadvance(S_new,U1,U2,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes2,fluxes3,THINCArr,parameters,dx,prob_lo,dt,time);

    FillPatch(*this, U.data, TWOGHOST, time, Phi_Type, 0, U1.data.nComp());

    U1 = ((U*(1.0/3.0))+(U2*(2.0/3.0)));


    MultiFab::Copy(S_new, U1.data, 0, 0, U1.data.nComp(), S_new.nGrow());
    FillPatch(*this, U1.data, TWOGHOST, time, Phi_Type, 0, U1.data.nComp());


    /*****************************************
     * This flux addition needs to be updated
     * to match the RK3 method
     * ***************************************/

    /*if(do_reflux)
    {
        for (int j = 0; j < BL_SPACEDIM; j++)
        {
            MultiFab::LinComb(fluxes[j],0.5,fluxes1[j],0,0.5,fluxes2[j],0,0,fluxes[j].nComp(),0);
        }
        ScaleFluxes(fluxes,dt,dx,accessPattern,parameters);
    }*/

}

void AmrLevelAdv::FourStage_RK3_HLLC_advance(MultiFab& S_new, Array<MultiFab, AMREX_SPACEDIM>& fluxes, ParameterStruct& parameters, const Real* dx, const Real* prob_lo, Real dt, Real time)
{

    const int TWOGHOST = 2;
    const int ONEGHOST = 1;

    /*-------------------------------------------------------------
   * Declared flux multifabs. Several are needed as we need to
   * keep track of the flux in each Runge-Kutta timestep to such
   * that we can reflux appropriately
   * -----------------------------------------------------------*/

    Array <MultiFab, AMREX_SPACEDIM> fluxes1;
    Array <MultiFab, AMREX_SPACEDIM> fluxes2;
    Array <MultiFab, AMREX_SPACEDIM> fluxes3;
    Array <MultiFab, AMREX_SPACEDIM> fluxes4;

    for (int j = 0; j < BL_SPACEDIM; j++)
    {
        BoxArray ba = S_new.boxArray();
        ba.surroundingNodes(j);

        fluxes1[j].define(ba, dmap, parameters.Ncomp, 0);
        fluxes2[j].define(ba, dmap, parameters.Ncomp, 0);
        fluxes3[j].define(ba, dmap, parameters.Ncomp, 0);
        fluxes4[j].define(ba, dmap, parameters.Ncomp, 0);

    }


    MultiFab S0(grids, dmap, parameters.Ncomp, TWOGHOST);
    FillPatch(*this, S0, TWOGHOST, time, Phi_Type, 0, parameters.Ncomp);

    MultiFab S1(grids, dmap, parameters.Ncomp, TWOGHOST);
    MultiFab S2(grids, dmap, parameters.Ncomp, TWOGHOST);
    MultiFab S3(grids, dmap, parameters.Ncomp, TWOGHOST);

    MultiFab SL(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SR(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SLStar(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SRStar(grids, dmap, parameters.Ncomp, ONEGHOST);
    MultiFab SStarStar(grids, dmap, parameters.Ncomp, ONEGHOST);

    THINCArray THINCArr(this->grids,this->dmap,ONEGHOST,parameters);

    CellArray U (S0,accessPattern,parameters);
    CellArray U1(S1,accessPattern,parameters);
    CellArray U2(S2,accessPattern,parameters);
    CellArray U3(S3,accessPattern,parameters);

    CellArray UL(SL,accessPattern,parameters);
    CellArray UR(SR,accessPattern,parameters);
    CellArray ULStar (SLStar,accessPattern,parameters);
    CellArray URStar (SRStar,accessPattern,parameters);
    CellArray UStarStar (SRStar,accessPattern,parameters);

    AMR_HLLCadvance(S_new,U,U2,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes2,fluxes3,THINCArr,parameters,dx,prob_lo,dt,time);

    U1 = ((U*(1.0/2.0))+(U2*(1.0/2.0)));

    AMR_HLLCadvance(S_new,U1,U3,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes2,fluxes3,THINCArr,parameters,dx,prob_lo,dt,time);

    U2 = ((U1*(1.0/2.0))+(U3*(1.0/2.0)));

    AMR_HLLCadvance(S_new,U2,U1,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes2,fluxes3,THINCArr,parameters,dx,prob_lo,dt,time);

    FillPatch(*this, U.data, TWOGHOST, time, Phi_Type, 0, U.data.nComp());

    U3 = ((U*(2.0/3.0))+(U2*(1.0/6.0))+(U1*(1.0/6.0)));

    AMR_HLLCadvance(S_new,U3,U2,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes2,fluxes3,THINCArr,parameters,dx,prob_lo,dt,time);

    U1 = ((U3*(1.0/2.0))+(U2*(1.0/2.0)));


    MultiFab::Copy(S_new, U1.data, 0, 0, U1.data.nComp(), S_new.nGrow());
    FillPatch(*this, U1.data, TWOGHOST, time, Phi_Type, 0, U1.data.nComp());


    /*****************************************
     * This flux addition needs to be updated
     * to match the RK3 method
     * ***************************************/

    /*if(do_reflux)
    {
        for (int j = 0; j < BL_SPACEDIM; j++)
        {
            MultiFab::LinComb(fluxes[j],0.5,fluxes1[j],0,0.5,fluxes2[j],0,0,fluxes[j].nComp(),0);
        }
        ScaleFluxes(fluxes,dt,dx,accessPattern,parameters);
    }*/

}

void AmrLevelAdv::AMR_HLLCadvance(MultiFab& S_new,CellArray& U,CellArray& U1, CellArray& UL, CellArray& UR, CellArray& ULStar, CellArray& URStar, CellArray& UStarStar, Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Array<MultiFab, AMREX_SPACEDIM>& flux_arr_diff, Array<MultiFab, AMREX_SPACEDIM>& flux_arr_Ustar,THINCArray& THINC,ParameterStruct& parameters, const Real* dx, const Real* prob_lo, Real dt, Real time)
{
    Direction_enum d;


    U.data.FillBoundary(geom.periodicity());
    FillDomainBoundary(U.data, geom, bc);

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

        if(parameters.MUSCL)
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for(MFIter mfi(UL.data); mfi.isValid(); ++mfi )
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

                MUSCLextrapolate(Ubox,ULbox,URbox,d);

            }


            if(parameters.THINC)
            {
#ifdef _OPENMP
#pragma omp parallel
#endif
                for(MFIter mfi(UL.data); mfi.isValid(); ++mfi )
                {
                    const Box& bx = mfi.validbox();

                    BoxAccessCellArray  Ubox(mfi,bx,U);
                    BoxAccessCellArray  ULbox(mfi,bx,UL);
                    BoxAccessCellArray  URbox(mfi,bx,UR);
                    BoxAccessCellArray ULTHINC(mfi,bx,ULStar);
                    BoxAccessCellArray URTHINC(mfi,bx,URStar);
                    BoxAccessTHINCArray THINCbox(mfi,bx,THINC);

                    THINCbox.THINCreconstruction(Ubox,ULbox,URbox,ULTHINC,URTHINC,parameters,dx,d);

                    ULbox.cleanUpAlpha(dx,prob_lo);
                    URbox.cleanUpAlpha(dx,prob_lo);

                }
            }


            UL.primitiveToConservative();
            UR.primitiveToConservative();

            UL.getSoundSpeed();
            UR.getSoundSpeed();

        }



        UL.data.FillBoundary(geom.periodicity());
        UR.data.FillBoundary(geom.periodicity());

        FillDomainBoundary(UL.data, geom, bc);
        FillDomainBoundary(UR.data, geom, bc);



        /*-------------------------------------------------------------
         * Calulate HLLC flux and update the new array.
         * -----------------------------------------------------------*/

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(S_new); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            FArrayBox& flux_fab1        = flux_arr[d][mfi];
            FArrayBox& flux_fab2        = flux_arr_diff[d][mfi];
            FArrayBox& flux_fab_UStar   = flux_arr_Ustar[d][mfi];

            BoxAccessCellArray Ubox(mfi,bx,U);
            BoxAccessCellArray U1box(mfi,bx,U1);
            BoxAccessCellArray ULbox(mfi,bx,UL);
            BoxAccessCellArray URbox(mfi,bx,UR);
            BoxAccessCellArray ULStarbox(mfi,bx,ULStar);
            BoxAccessCellArray URStarbox(mfi,bx,URStar);
            BoxAccessCellArray UStarStarbox(mfi,bx,UStarStar);
            BoxAccessCellArray fluxbox1(bx,flux_fab1,U);
            BoxAccessCellArray fluxbox2(bx,flux_fab2,U);
            BoxAccessCellArray fluxStyleUstar(bx,flux_fab_UStar,U);


            if(parameters.SOLID)
            {
                calc_5Wave_fluxes(fluxbox1, ULbox, URbox, ULStarbox, URStarbox,fluxStyleUstar, parameters,d,dx,prob_lo); //UStarStarbox
            }
            else
            {
                calc_fluxes(fluxbox1, ULbox, URbox, fluxStyleUstar, parameters,d,dx,prob_lo);
            }

            update(fluxbox1, Ubox, U1box, parameters,d,dt,dx);

            //calculatePathConservativeFluxes(fluxbox1,fluxbox2,Ubox, ULbox, URbox,fluxStyleUstar,parameters,d,dx);

            //PCupdate(fluxbox1,fluxbox2, Ubox, U1box, parameters,d,dt,dx);

            U1box.conservativeToPrimitive();

            if(parameters.REACTIVE)
            {
                reactiveUpdateInHLLC(U1box,parameters,dt);
            }

        }
    }



    U1.data.FillBoundary(geom.periodicity());
    FillDomainBoundary(U1.data, geom, bc);

    U1.conservativeToPrimitive();

}

void AmrLevelAdv::AMR_PathConservative_HLLCadvance(MultiFab& S_new,CellArray& U,CellArray& U1, CellArray& UL, CellArray& UR, CellArray& ULStar, CellArray& URStar, CellArray& UStarStar, Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Array<MultiFab, AMREX_SPACEDIM>& flux_arr_diff, Array<MultiFab, AMREX_SPACEDIM>& flux_arr_Ustar,THINCArray& THINC,ParameterStruct& parameters, const Real* dx, const Real* prob_lo, Real dt, Real time)
{
    Direction_enum d; 


    U.data.FillBoundary(geom.periodicity());
    FillDomainBoundary(U.data, geom, bc);

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

        if(parameters.MUSCL)
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for(MFIter mfi(UL.data); mfi.isValid(); ++mfi )
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

                MUSCLextrapolate(Ubox,ULbox,URbox,d);

            }


            if(parameters.THINC)
            {
#ifdef _OPENMP
#pragma omp parallel
#endif
                for(MFIter mfi(UL.data); mfi.isValid(); ++mfi )
                {
                    const Box& bx = mfi.validbox();

                    BoxAccessCellArray  Ubox(mfi,bx,U);
                    BoxAccessCellArray  ULbox(mfi,bx,UL);
                    BoxAccessCellArray  URbox(mfi,bx,UR);
                    BoxAccessCellArray ULTHINC(mfi,bx,ULStar);
                    BoxAccessCellArray URTHINC(mfi,bx,URStar);
                    BoxAccessTHINCArray THINCbox(mfi,bx,THINC);

                    THINCbox.THINCreconstruction(Ubox,ULbox,URbox,ULTHINC,URTHINC,parameters,dx,d);

                    ULbox.cleanUpAlpha(dx,prob_lo);
                    URbox.cleanUpAlpha(dx,prob_lo);

                }
            }


            UL.primitiveToConservative();
            UR.primitiveToConservative();

            UL.getSoundSpeed();
            UR.getSoundSpeed();

        }



        UL.data.FillBoundary(geom.periodicity());
        UR.data.FillBoundary(geom.periodicity());

        FillDomainBoundary(UL.data, geom, bc);
        FillDomainBoundary(UR.data, geom, bc);

        

        /*-------------------------------------------------------------
         * Calulate HLLC flux and update the new array.
         * -----------------------------------------------------------*/

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(S_new); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            FArrayBox& flux_fab1        = flux_arr[d][mfi];
            FArrayBox& flux_fab2        = flux_arr_diff[d][mfi];
            FArrayBox& flux_fab_UStar   = flux_arr_Ustar[d][mfi];

            BoxAccessCellArray Ubox(mfi,bx,U);
            BoxAccessCellArray U1box(mfi,bx,U1);
            BoxAccessCellArray ULbox(mfi,bx,UL);
            BoxAccessCellArray URbox(mfi,bx,UR);
            BoxAccessCellArray ULStarbox(mfi,bx,ULStar);
            BoxAccessCellArray URStarbox(mfi,bx,URStar);
            BoxAccessCellArray UStarStarbox(mfi,bx,UStarStar);
            BoxAccessCellArray fluxbox1(bx,flux_fab1,U);
            BoxAccessCellArray fluxbox2(bx,flux_fab2,U);
            BoxAccessCellArray fluxStyleUstar(bx,flux_fab_UStar,U);


            if(parameters.SOLID)
            {
                calc_5Wave_fluxes(fluxbox1, ULbox, URbox, ULStarbox, URStarbox, fluxStyleUstar, parameters,d,dx,prob_lo); //UStarStarbox
            }
            else
            {
                calc_fluxes(fluxbox1, ULbox, URbox, fluxStyleUstar, parameters,d,dx,prob_lo);
            }

            update(fluxbox1, Ubox, U1box, parameters,d,dt,dx);

            //calculatePathConservativeFluxes(fluxbox1,fluxbox2,Ubox, ULbox, URbox,fluxStyleUstar,parameters,d,dx);

            //PCupdate(fluxbox1,fluxbox2, Ubox, U1box, parameters,d,dt,dx);

            //U1box.conservativeToPrimitive();

        }
    }

    U1.conservativeToPrimitive();


    if(parameters.REACTIVE)
    {
        for(MFIter mfi(S_new); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();
            BoxAccessCellArray U1box(mfi,bx,U1);

            reactiveUpdateInHLLC(U1box,parameters,dt);
        }
    }



    U1.data.FillBoundary(geom.periodicity());
    FillDomainBoundary(U1.data, geom, bc);

    U1.conservativeToPrimitive();

}

Real AmrLevelAdv::advance (Real time, Real dt, int  iteration, int  ncycle)
{

    /*-------------------------------------------------------------
     * Runge-Kutta time integration with MUSCL/HLLC is used to update the
     * solution on a given level.
     * -----------------------------------------------------------*/

    MultiFab& S_new = get_new_data(Phi_Type);



    const Real prev_time = state[Phi_Type].prevTime();
    const Real cur_time = state[Phi_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();



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

    /*-------------------------------------------------------------
     * Declared flux multifabs. Several are needed as we need to
     * keep track of the flux in each Runge-Kutta timestep to such
     * that we can reflux appropriately
     * -----------------------------------------------------------*/

    Array <MultiFab, AMREX_SPACEDIM> fluxes;
    Array <MultiFab, AMREX_SPACEDIM> fluxes1;
    Array <MultiFab, AMREX_SPACEDIM> fluxes2;
    Array <MultiFab, AMREX_SPACEDIM> UStarflux;

    for (int j = 0; j < AMREX_SPACEDIM; j++)
    {
        BoxArray ba = S_new.boxArray();
        ba.surroundingNodes(j);
        fluxes [j].define(ba, dmap, parameters.Ncomp, 0);
        fluxes1[j].define(ba, dmap, parameters.Ncomp, 0);
        fluxes2[j].define(ba, dmap, parameters.Ncomp, 0);
        UStarflux[j].define(ba, dmap, parameters.Ncomp, 0);
    }


    /*-------------------------------------------------------------
     * Declare some multifabs with Ghost cells to hold intermediate
     * data in the RK update, and then wrap them in CellArrays
     * -----------------------------------------------------------*/

    /*------------------------------------------------------------
     * States with Two Ghost Cells:
     * -----------------------------------------------------------*/
    int TWOGHOST = 2;

     MultiFab S0(grids, dmap, parameters.Ncomp, TWOGHOST);
     FillPatch(*this, S0, TWOGHOST, time, Phi_Type, 0, parameters.Ncomp);
     MultiFab S1(grids, dmap, parameters.Ncomp, TWOGHOST);
     FillPatch(*this, S1, TWOGHOST, time, Phi_Type, 0, parameters.Ncomp);
     MultiFab S2(grids, dmap, parameters.Ncomp, TWOGHOST);
     FillPatch(*this, S2, TWOGHOST, time, Phi_Type, 0, parameters.Ncomp);

     /*-------------------------------------------------------------
      * States with One Ghost Cell:
      * -----------------------------------------------------------*/

     int ONEGHOST = 1;

     MultiFab SL(grids, dmap, parameters.Ncomp, ONEGHOST);
     FillPatch(*this, SL, ONEGHOST, time, Phi_Type, 0, parameters.Ncomp);
     MultiFab SR(grids, dmap, parameters.Ncomp, ONEGHOST);
     FillPatch(*this, SR, ONEGHOST, time, Phi_Type, 0, parameters.Ncomp);
     MultiFab SLStar(grids, dmap, parameters.Ncomp, ONEGHOST);
     FillPatch(*this, SLStar, ONEGHOST, time, Phi_Type, 0, parameters.Ncomp);
     MultiFab SRStar(grids, dmap, parameters.Ncomp, ONEGHOST);
     FillPatch(*this, SRStar, ONEGHOST, time, Phi_Type, 0, parameters.Ncomp);
     MultiFab SStarStar(grids, dmap, parameters.Ncomp, ONEGHOST);
     FillPatch(*this, SStarStar, ONEGHOST, time, Phi_Type, 0, parameters.Ncomp);
     MultiFab Sgrad(grids, dmap, parameters.Ncomp, ONEGHOST);
     FillPatch(*this, Sgrad, ONEGHOST, time, Phi_Type, 0, parameters.Ncomp);

    CellArray U (S0,accessPattern,parameters);
    CellArray U1(S1,accessPattern,parameters);
    CellArray U2(S2,accessPattern,parameters);
    CellArray U3(S2,accessPattern,parameters);
    CellArray UL(SL,accessPattern,parameters);
    CellArray UR(SR,accessPattern,parameters);
    CellArray ULStar(SLStar,accessPattern,parameters);
    CellArray URStar(SRStar,accessPattern,parameters);
    CellArray UStarStar(SStarStar,accessPattern,parameters);

    THINCArray THINCArr(this->grids,this->dmap,ONEGHOST,parameters);


    if(parameters.RADIAL)
    {
        geometricSourceTerm(U,parameters,dx,dt/2.0,prob_lo,S_new);
    }

    AMR_HLLCadvance(S_new,U ,U1,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes,UStarflux,THINCArr,parameters,dx,prob_lo,dt,time);

    AMR_HLLCadvance(S_new,U1,U2,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes,UStarflux,THINCArr,parameters,dx,prob_lo,dt,time);

    U1 = ((U*(3.0/4.0))+(U2*(1.0/4.0)));

    AMR_HLLCadvance(S_new,U1,U2,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes,UStarflux,THINCArr,parameters,dx,prob_lo,dt,time);

    FillPatch(*this, U.data, TWOGHOST, time, Phi_Type, 0, U1.data.nComp());

    U1 = ((U*(1.0/3.0))+(U2*(2.0/3.0)));


    /*AMR_HLLCadvance(S_new,U ,U1,UL,UR,ULStar,URStar,UStarStar,fluxes1,fluxes,UStarflux,THINCArr,parameters,dx,prob_lo,dt,time);

    AMR_HLLCadvance(S_new,U1,U2,UL,UR,ULStar,URStar,UStarStar,fluxes2,fluxes,UStarflux,THINCArr,parameters,dx,prob_lo,dt,time);

    U1 = ((U*(1.0/2.0))+(U2*(1.0/2.0)));*/

    /*if(parameters.REACTIVE)
    {
        reactiveUpdateOutHLLC(U1,U,U2,parameters,dt);
    }*/

    if(parameters.PLASTIC)
    {
        plastic.plasticUpdate(U1,parameters,dt,S_new);
    }

    if(parameters.RADIAL)
    {
        geometricSourceTerm(U1,parameters,dx,dt/2.0,prob_lo,S_new);
    }

    if(parameters.SOLID)
    {
        U1.cleanUpV();
    }

    {
        U1.cleanUpAlpha(dx,prob_lo);
    }


    U1.data.FillBoundary(geom.periodicity());
    FillDomainBoundary(U1.data, geom, bc);

    MultiFab::Copy(S_new, U1.data, 0, 0, U1.data.nComp(), S_new.nGrow());
    FillPatch(*this, U1.data, TWOGHOST, time, Phi_Type, 0, parameters.Ncomp);

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


    return dt;
}


Real AmrLevelAdv::new_advance (Real time, Real dt, int  iteration, int  ncycle)
{

    /*-------------------------------------------------------------
     * Runge-Kutta time integration with MUSCL/HLLC is used to update the
     * solution on a given level.
     * -----------------------------------------------------------*/

    MultiFab& S_new = get_new_data(Phi_Type);

    CellArray U1(S_new,accessPattern,parameters);

    const Real prev_time = state[Phi_Type].prevTime();
    const Real cur_time = state[Phi_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();


    /*********************************************************************
    * Get pointers to Flux registers, or set pointer to zero if not there.
    *********************************************************************/

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

    /**************************************
     * Array to store fluxes before adding
     * them to the flux register
     *************************************/

    Array <MultiFab, AMREX_SPACEDIM> fluxes ;

    for (int j = 0; j < BL_SPACEDIM; j++)
    {
        BoxArray ba = S_new.boxArray();

        ba.surroundingNodes(j);

        fluxes[j].define(ba, dmap, parameters.Ncomp, 0);
    }


    if(parameters.RADIAL)
    {
        geometricSourceTerm(U1,parameters,dx,dt/2.0,prob_lo,S_new);
    }


    /******************************
     * Choose update order
     *****************************/

    //SINGLE_HLLC_advance(S_new,fluxes,parameters,dx,prob_lo,dt,time);

    RK2_HLLC_advance(S_new,fluxes,parameters,dx,prob_lo,dt,time);

    //RK3_HLLC_advance(S_new,fluxes,parameters,dx,prob_lo,dt,time);

    //FourStage_RK3_HLLC_advance(S_new,fluxes,parameters,dx,prob_lo,dt,time);

    /******************************/

    //sourceUpdate(dt,mat);

    /*if(parameters.REACTIVE)
    {
        #ifdef _OPENMP
        #pragma omp parallel
        #endif
        for(MFIter mfi(S_new); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.validbox();

            BoxAccessCellArray U1box(mfi,bx,U1);

            reactiveUpdateInHLLC(U1box,parameters,dt);

        }

        U1.conservativeToPrimitive();
    }*/


    /*if(parameters.RADIAL)
    {
        geometricSourceTerm(U,parameters,dx,dt/2.0,prob_lo,S_new);
    }*/
    

    /*if(parameters.REACTIVE)
    {
        reactiveUpdate(U,U1,U2,parameters,dt,S_new);
    }*/

    if(parameters.PLASTIC)
    {
        plastic.plasticUpdate(U1,parameters,dt,S_new);
    }

    if(parameters.RADIAL)
    {
        geometricSourceTerm(U1,parameters,dx,dt/2.0,prob_lo,S_new);
    }

    if(parameters.SOLID)
    {
        U1.cleanUpV();
    }


    {
        U1.cleanUpAlpha(dx,prob_lo);
    }


    U1.data.FillBoundary(geom.periodicity());
    FillDomainBoundary(U1.data, geom, bc);
    FillPatch(*this, U1.data, 0, time, Phi_Type, 0, parameters.Ncomp);

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

    return dt;
}

Real AmrLevelAdv::estTimeStep (Real)
{

    /*-------------------------------------------------------------
     * This timestep function is essentially unchanged from the
     * tutorial execept that we calculate the wavespeed in the
     * normal way, instead of the plain advection way.
     * -----------------------------------------------------------*/

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

    if(nStep() <= 5)
    {
        dt_est *= 0.2;
    }

    if (verbose) {
	amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level 
                       << ":  dt_est = " << dt_est << std::endl;
    }
    
    return dt_est;
}

void findBiggestGradientOnLevel(const Box& box, BoxAccessCellArray& U, BoxAccessCellArray& grad,const Real* dx, Real time, Real level, Real& gradMax, MaterialSpecifier& n)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Real ax,ay;

    for 		(int k = lo.z; k <= hi.z; ++k)
    {
        for 	(int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                ax = (U(i+1,j,k,n)-U(i-1,j,k,n))/(2.0*dx[0]);
                ay = (U(i,j+1,k,n)-U(i,j-1,k,n))/(2.0*dx[1]);

                grad(i,j,k,n) = std::max(std::abs(ax),std::abs(ay)); //sqrt(ax*ax+ay*ay);

                gradMax = (grad(i,j,k,n) > gradMax ? grad(i,j,k,n) : gradMax);
            }
        }
    }

    return;
}

void findBiggestDifferenceOnLevel(const Box& box, BoxAccessCellArray& U, BoxAccessCellArray& diff,const Real* dx, Real time, Real level, Real& diffMax, MaterialSpecifier& n)
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real ax,ay;

    for 		(int k = lo.z; k <= hi.z; ++k)
    {
        for 	(int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                ax = (U(i+1,j,k,n)-U(i-1,j,k,n));
                ay = (U(i,j+1,k,n)-U(i,j-1,k,n));

                diff(i,j,k,n) = (std::max(ax*ax,ay*ay))/(U(i,j,k,n)*U(i,j,k,n)); //sqrt(ax*ax+ay*ay);

                diffMax = (diff(i,j,k,n) > diffMax ? diff(i,j,k,n) : diffMax);
            }
        }
    }

    return;
}

void C_state_error_grad(Array4<char> const& tagarr, const Box& box, BoxAccessCellArray& U, BoxAccessCellArray& grad, const Real* dx, Real time, Real level, Vector<Real>& gradMax, Vector<Real> levelGradientCoefficients, AccessPattern& accessPattern)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);

    int m = 0;

    for(auto n : accessPattern.refineVariables)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    if(grad(i,j,k,n) > levelGradientCoefficients[level]*gradMax[m])
                    {
                        tagarr(i,j,k) = TagBox::SET;
                    }
                }
            }
        }

        m++;
    }

	return;

}

void C_state_error_diff(Array4<char> const& tagarr, const Box& box, BoxAccessCellArray& U, BoxAccessCellArray& grad, const Real* dx, Real time, Real level, Vector<Real>& gradMax, Vector<Real> levelGradientCoefficients, AccessPattern& accessPattern)
{
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real tolerance = 1E-10;

    int m = 0;

    for(auto n : accessPattern.refineVariables)
    {
        for 		(int k = lo.z; k <= hi.z; ++k)
        {
            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    if(grad(i,j,k,n) > levelGradientCoefficients[level] && grad(i,j,k,n) > tolerance)
                    {
                        tagarr(i,j,k) = TagBox::SET;
                    }
                }
            }
        }

        m++;
    }

    return;

}

void AmrLevelAdv::errorEst (TagBoxArray& tags, int clearval, int tagval, Real time, int n_error_buf, int ngrow)
{

    /*-------------------------------------------------------------
     * This function tags cells for regridding. We first calculate
     * the maximum gradient on a level, and tag cells based on
     * whether they are within a certain fraction of that gradient
     * -----------------------------------------------------------*/

    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(Phi_Type);
    
    // State with ghost cells
    MultiFab Sborder(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);

    MultiFab Sgrad(grids, dmap, parameters.Ncomp, NUM_GROW);
    FillPatch(*this, Sgrad, NUM_GROW, time, Phi_Type, 0, parameters.Ncomp);
    
    CellArray U(Sborder,accessPattern,parameters);
    CellArray grad(Sgrad,accessPattern,parameters);

    Vector<Real> gradMaxVec(parameters.Ncomp,0.0); //accessPattern.refineVariables.size()+1

    int m = 0;

    for(auto n : accessPattern.refineVariables)
    {
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
                BoxAccessCellArray 			gradbox(mfi,bx,grad);

                findBiggestDifferenceOnLevel(bx,Ubox,gradbox,dx,time,level,gradMax,n);
                //findBiggestGradientOnLevel(bx,Ubox,gradbox,dx,time,level,gradMax,n);
            }
        }

        ParallelDescriptor::ReduceRealMax(gradMax);

        gradMaxVec[m] = gradMax;

        m++;
    }
	

	
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
		for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
		{
			const Box& bx	  			= mfi.validbox();
			
            BoxAccessCellArray 			Ubox(mfi,bx,U);
            BoxAccessCellArray 			gradbox(mfi,bx,grad);

            TagBox&    tagfab	  		= tags[mfi];
            
            Array4<char> const& tagarr 	= tagfab.array();
            
            C_state_error_diff(tagarr,bx,Ubox,gradbox,dx,time,level,gradMaxVec,levelGradientCoefficients,accessPattern);
            //C_state_error_grad(tagarr,bx,Ubox,gradbox,dx,time,level,gradMaxVec,levelGradientCoefficients,accessPattern);

        }
	}
}

void AmrLevelAdv::variableSetUp ()
{

    BL_ASSERT(desc_lst.size() == 0);

    read_params();

    desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),StateDescriptor::Point,0,parameters.Ncomp,&cell_cons_interp);

    bc.resize(parameters.Ncomp);

    setBoundaryConditions(bc,parameters,initial,accessPattern);

    for(int n = 0; n<parameters.Ncomp; n++)
    {
        desc_lst.setComponent(Phi_Type, n, accessPattern.variableNames[n] , bc[n], StateDescriptor::BndryFunc(phifill));
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
    pp.queryarr("level_Coeff",levelGradientCoefficients);
    
    libConfigInitialiseDataStructs(parameters,initial,plastic);
    
    accessPattern.define(parameters);

    /*for(auto n : accessPattern.refineVariables)
    {
        Print() << accessPattern.variableNames[accessPattern[n.var]+n.mat] << std::endl;
    }*/
    
    parameters.Ncomp = accessPattern.variableNames.size();

    Geometry const* gg = AMReX::top()->getDefaultGeometry();

    // This tutorial code only supports Cartesian coordinates.
    if (! gg->IsCartesian()) {
	amrex::Abort("Please set geom.coord_sys = 0");
    }
}

#include "simulationheader.h"

using namespace amrex;


void main_main ()
{
    Real strt_time = amrex::second();

    int max_grid_size, nsteps_max, plot_int;

    Array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(0,0,0)};  // non-periodic in all directions by default

    ParmParse pp;

    ParameterStruct parameters;
    InitialStruct   initial;

    initialiseDataStructs(parameters,initial);

    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(parameters.n_cells[0]-1, parameters.n_cells[1]-1, parameters.n_cells[2]-1));

    Box domain(dom_lo, dom_hi);

    BoxArray ba(domain);

    ba.maxSize(parameters.max_grid_size);

    RealBox real_box({AMREX_D_DECL(0.0,0.0,0.0)},{AMREX_D_DECL(parameters.dimL[0], parameters.dimL[1], parameters.dimL[2])});

    Geometry geom(domain,real_box,CoordSys::cartesian,is_periodic);

    DistributionMapping dm(ba);

    MultiFab phi_old(ba, dm, parameters.Ncomp, parameters.Nghost);
    MultiFab phi_new(ba, dm, parameters.Ncomp, parameters.Nghost);
    MultiFab phi_Star(ba, dm, parameters.Ncomp, parameters.Nghost);

    TimeStep timeStep(ba, dm, AMREX_SPACEDIM, parameters.Nghost);

    /*  Components:
     * 0 rho
     * 1 rhou
     * 2 E
     * 3 u
     * 4 p
     */


    Vector<BCRec> bc(phi_new.nComp());

    for(int n = 0; n < phi_new.nComp(); ++n)
    {
        for(int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bc[n].setLo(idim, BCType::foextrap);
            bc[n].setHi(idim, BCType::foextrap);
        }
    }


    setInitialConditions(phi_new,parameters,initial);


    /*

    if (parameters.plot_int > 0)
    {
        int n = 0;

        //WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, 0);

        std::string pltfile = amrex::Concatenate("data/rho",n,1);
        WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, 0.0, 0);
        PrintTo1DGnuplotFile(phi_new, pltfile, RHO);
        pltfile = amrex::Concatenate("data/u",n,1);
        PrintTo1DGnuplotFile(phi_new, pltfile, VELOCITY);
        pltfile = amrex::Concatenate("data/p",n,1);
        PrintTo1DGnuplotFile(phi_new, pltfile, P);
    }


    Array<MultiFab, AMREX_SPACEDIM> flux_arr;

    for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        BoxArray edge_ba = ba;

        edge_ba.surroundingNodes(dir);

        flux_arr[dir].define(edge_ba, dm, parameters.Ncomp, 0);
    }

    parameters.dx[0] = geom.CellSize()[0];
    parameters.dx[1] = geom.CellSize()[1];
    parameters.dx[2] = geom.CellSize()[2];

    int n = 0;


    for(Real t = 0.0 ; t<initial.finalT; t += parameters.dt, n++)
    {

        FillDomainBoundary(phi_new, geom, bc);
        phi_new.FillBoundary(geom.periodicity());

        parameters.dt = getTimeStep(phi_new,timeStep,parameters);

        //amrex::Print() << "dt: " << parameters.dt << "      t: " << t << " / " << initial.finalT << std::endl;

        if(t+parameters.dt>initial.finalT)
        {
            parameters.dt=initial.finalT-t;
        }



        MultiFab::Copy(phi_old, phi_new, 0, 0, parameters.Ncomp, parameters.Nghost);

        advance(phi_old, phi_new, phi_Star, flux_arr, geom, parameters);

    }

    std::string pltfile = amrex::Concatenate("data/rho",1,1);
    WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, initial.finalT, 0);
    PrintTo1DGnuplotFile(phi_new, pltfile, RHO);
    pltfile = amrex::Concatenate("data/u",1,1);
    PrintTo1DGnuplotFile(phi_new, pltfile, VELOCITY);
    pltfile = amrex::Concatenate("data/p",1,1);
    PrintTo1DGnuplotFile(phi_new, pltfile, P);

    pltfile = amrex::Concatenate("data/rhoU",1,1);
    PrintTo1DGnuplotFile(phi_new, pltfile, RHOU);
    pltfile = amrex::Concatenate("data/E",1,1);
    PrintTo1DGnuplotFile(phi_new, pltfile, TOTAL_E);

    */

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = amrex::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();

    return 0;
}

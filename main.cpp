#include "simulationheader.h"
#include "timestep.h"
using namespace amrex;


void main_main ()
{
    Real start_time = amrex::second();

    int max_grid_size, nsteps_max, plot_int;

    Array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(0,0,0)};  // non-periodic in all directions by default

    ParmParse pp;

    ParameterStruct parameters;
    InitialStruct   initial;

    initialiseDataStructs(parameters,initial);

    Print() << initial.filename << std::endl;

    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(parameters.n_cells[0]-1, parameters.n_cells[1]-1, parameters.n_cells[2]-1));

    Box domain(dom_lo, dom_hi);

    BoxArray ba(domain);

    ba.maxSize(parameters.max_grid_size);

    RealBox real_box({AMREX_D_DECL(0.0,0.0,0.0)},{AMREX_D_DECL(parameters.dimL[0], parameters.dimL[1], parameters.dimL[2])});

    Geometry geom(domain,real_box,CoordSys::cartesian,is_periodic);

    DistributionMapping dm(ba);

    CellArray U(ba, dm, parameters.Ncomp, parameters.Nghost);
    CellArray U1(ba, dm, parameters.Ncomp, parameters.Nghost);
    CellArray U2(ba, dm, parameters.Ncomp, parameters.Nghost);
    CellArray UL(ba, dm, parameters.Ncomp, parameters.Nghost);
    CellArray UR(ba, dm, parameters.Ncomp, parameters.Nghost);
    CellArray UStar(ba, dm, parameters.Ncomp, parameters.Nghost);
    CellArray MUSCLgrad(ba, dm, parameters.Ncomp, parameters.Nghost);

    TimeStep timeStep(ba, dm, AMREX_SPACEDIM, parameters.Nghost);

    Vector<BCRec> bc(U1.data.nComp());

    for(int n = 0; n < U1.data.nComp(); ++n)
    {
        for(int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bc[n].setLo(idim, BCType::foextrap);
            bc[n].setHi(idim, BCType::foextrap);
        }
    }

    setInitialConditions(U1,parameters);

    PrintAllVarsTo1DGnuplotFile(U1.data,0,initial.filename);

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

        FillDomainBoundary(U1.data, geom, bc);

        U1.data.FillBoundary(geom.periodicity());

        parameters.dt = timeStep.getTimeStep(U1,parameters);

        amrex::Print() << "dt: " << parameters.dt << "      t: " << t << " / " << initial.finalT << std::endl;

        if(t+parameters.dt>initial.finalT)
        {
            parameters.dt=initial.finalT-t;
        }

        U = U1;

        RKadvance(U, U1, U2, UL, UR, MUSCLgrad, UStar, flux_arr, geom, parameters,bc);

    }

    PrintAllVarsTo1DGnuplotFile(U1.data,1,initial.filename);



    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = amrex::second() - start_time;
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

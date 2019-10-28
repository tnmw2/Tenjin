#include "simulationheader.h"
#include "timestep.h"
using namespace amrex;

void main_main ()
{
    Real start_time = amrex::second();

    ParmParse pp;

    /* -----------------------------------------------------
     * parameters and initial are structs that hold and pass
     * around simulation parameters.
     * -----------------------------------------------------*/

    ParameterStruct parameters;
    InitialStruct   initial;
    PlasticEOS      plastic;

    libConfigInitialiseDataStructs(parameters,initial,plastic);

    /* ----------------------------------------------------
     * Declare the domain and geometry required for Multifabs
     * Also delcare a Distribution mapping (Only the standard
     * one has been used, but I would like to know the benefits
     * of using others)
     * -----------------------------------------------------*/

    IntVect lo(AMREX_D_DECL(0,0,0));
    IntVect hi(AMREX_D_DECL(parameters.n_cells[0]-1, parameters.n_cells[1]-1, parameters.n_cells[2]-1));

    Box domain(lo,hi);

    BoxArray ba(domain);

    ba.maxSize(parameters.max_grid_size);

    RealBox real_box({AMREX_D_DECL(0.0,0.0,0.0)},{AMREX_D_DECL(parameters.dimL[0], parameters.dimL[1], parameters.dimL[2])});

    Geometry geom(domain,real_box,CoordSys::cartesian,{AMREX_D_DECL(0,0,0)});

    DistributionMapping dm(ba);

    parameters.dx = {geom.CellSize()[0],geom.CellSize()[1],geom.CellSize()[2]};

    /* ----------------------------------------------------
     *  Because the number of components in a Multifab changes
     * from simulation to simulation, we declare an
     * AccessPattern which tells us how to get a specific
     * variable from a Multifab.
     * ------------------------------------------------------*/

    AccessPattern accessPattern(parameters);

    parameters.Ncomp = accessPattern.variableNames.size();

    /* ----------------------------------------------------
     * Declare Multifab wrappers with overloaded =,+,* for
     * ease.
     * Declare a timestep class that will calculate the
     * timestep and have its own Multifab to store data.
     * ----------------------------------------------------*/

    CellArray U         (ba, dm, parameters.Ncomp, parameters.Nghost, accessPattern, parameters);
    CellArray U1        (ba, dm, parameters.Ncomp, parameters.Nghost, accessPattern, parameters);
    CellArray U2        (ba, dm, parameters.Ncomp, parameters.Nghost, accessPattern, parameters);
    CellArray UL        (ba, dm, parameters.Ncomp, parameters.Nghost, accessPattern, parameters);
    CellArray UR        (ba, dm, parameters.Ncomp, parameters.Nghost, accessPattern, parameters);
    CellArray ULStar    (ba, dm, parameters.Ncomp, parameters.Nghost, accessPattern, parameters);
    CellArray URStar    (ba, dm, parameters.Ncomp, parameters.Nghost, accessPattern, parameters);
    CellArray UStarStar (ba, dm, parameters.Ncomp, parameters.Nghost, accessPattern, parameters);
    CellArray MUSCLgrad (ba, dm, parameters.Ncomp, parameters.Nghost, accessPattern, parameters);

    THINCArray THINCArr (ba, dm,                   parameters.Nghost,                parameters);

    TimeStep timeStep(ba, dm, AMREX_SPACEDIM, parameters.Nghost);

    /* ----------------------------------------------------
     * Declare transmissive boundary conditions and set up
     * a flux vector
     * ----------------------------------------------------*/

    Array<MultiFab, AMREX_SPACEDIM> flux_arr;

    for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        BoxArray edge_ba = ba;

        edge_ba.surroundingNodes(dir);

        flux_arr[dir].define(edge_ba, dm, parameters.Ncomp, 0);
    }

        Vector<BCRec> bc(parameters.Ncomp);
        setBoundaryConditions(bc,parameters,initial,accessPattern);

    /* ----------------------------------------------------
     * Setup and Print the initial conditions
     * ----------------------------------------------------*/

    setInitialConditions(U1,parameters,initial);

    {
        const std::string& pltfile = Concatenate(initial.filename,0,5);

        WriteSingleLevelPlotfile(pltfile, U1.data, U1.accessPattern.variableNames , geom, 0.0, 0);
    }


    /* ----------------------------------------------------
     * The main time loop
     * ----------------------------------------------------*/


    //return;

    int n = 0;

    int take_pic_counter = 0;

    for(Real t = 0.0 ; t<initial.finalT; t += parameters.dt, n++)
    {

        FillDomainBoundary(U1.data, geom, bc);

        U1.data.FillBoundary(geom.periodicity());

        parameters.dt = timeStep.getTimeStep(U1,parameters);

        if(n%1==0)
        {
            amrex::Print() << "dt: " << parameters.dt << "      t: " << t << " / " << initial.finalT << std::endl;
        }

        if(t+parameters.dt>initial.finalT)
        {
            parameters.dt=initial.finalT-t;
        }

        if(n<5)
        {
            parameters.dt*=0.2;
        }

        U = U1;

        advance(U, U1, U2, UL, UR, MUSCLgrad, ULStar, URStar, UStarStar, flux_arr, geom, parameters,bc,THINCArr);

        if(parameters.REACTIVE)
        {
            reactiveUpdate(U,U1,U2,parameters);
        }

        if(parameters.PLASTIC)
        {
            plastic.plasticUpdate(U1,parameters);
        }

        if(parameters.RADIAL)
        {
            geometricSourceTerm(U1,parameters);
        }

        if(parameters.SOLID)
        {
            U1.cleanUpV();
        }

        if(t > ((Real)take_pic_counter)*(initial.finalT)/((Real)initial.numberOfPictures))
        {
            take_pic_counter++;
            const std::string& pltfile = Concatenate(initial.filename,n,5);

            WriteSingleLevelPlotfile(pltfile, U1.data, U1.accessPattern.variableNames , geom, t, n);

        }


        if(U1.data.contains_nan())
        {
            Print() << "Nan found in U1" << std::endl;
            break;
        }

        //break;

    }



    {
        const std::string& pltfile = Concatenate(initial.filename,n,5);
        //PrintAllVarsTo1DGnuplotFile(U1,1,initial.filename);
        WriteSingleLevelPlotfile(pltfile, U1.data, U1.accessPattern.variableNames , geom, initial.finalT, n);
    }


    Real stop_time = amrex::second() - start_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);
    Print() << "Run time = " << stop_time << std::endl;

}

int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();

    return 0;
}

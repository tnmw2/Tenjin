#include <new>
#include <iostream>
#include <iomanip>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>

#include "simulationheader.h"

void stringtochar(std::string string,char* file);


using namespace amrex;

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    Real dRunTime1 = amrex::second();

    int  max_step;
    Real strt_time;
    Real stop_time;

    {
        using namespace libconfig;

        ParmParse pp;

        max_step  = -1;
        strt_time =  0.0;
        stop_time = -1.0;

        pp.query("max_step",max_step);
        pp.query("strt_time",strt_time);
        pp.query("stop_time",stop_time);

        std::string settingsFileString;

        pp.get("SettingsFile",settingsFileString);

        char settingsFile[settingsFileString.length()+1];

        stringtochar(settingsFileString,settingsFile);

        Config cfg;

        try
        {
            cfg.readFile(settingsFile);
        }
        catch(ParseException except)
        {
            std::cout << "Incorrect Settings file" << std::endl;
            exit(1);
        }

        cfg.lookupValue("finalTime",stop_time);
    }

    if (strt_time < 0.0) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0) {
    amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {

    Amr amr;

    amr.init(strt_time,stop_time);

    int a = 0;

    while ( amr.okToContinue() &&
           (amr.levelSteps(0) < max_step || max_step < 0) &&
           (amr.cumTime() < stop_time || stop_time < 0.0) )

    {
        //
        // Do a coarse timestep.  Recursively calls timeStep()
        //
        amr.coarseTimeStep(stop_time);



        a++;
    }

    // Write final checkpoint and plotfile
    if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
        amr.checkPoint();
    }

    if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
        amr.writePlotFile();
    }

    }

    Real dRunTime2 = amrex::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;

    amrex::Finalize();

    return 0;
}

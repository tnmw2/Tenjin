#include "simulationheader.h"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    //main_main();

    amrex::Finalize();

    std::cout << "Something is working" << std::endl;

    return 0;
}

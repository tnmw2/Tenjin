#include "simulationheader.h"

void PrintTo1DGnuplotFile(CellArray& U, std::string filename, std::string name, int picture, Variable var, int mat =0)
{
    std::string pltfilevarstub = filename;

    std::string pltfilevar = pltfilevarstub.append(name);

    const std::string pltfilefinal = amrex::Concatenate(pltfilevar,picture,2);

    for( MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& box = mfi.validbox();

        BoxAccessCellArray Ubox(mfi,box,U);

        const auto lo = lbound(Ubox.box);
        const auto hi = ubound(Ubox.box);

        for    		(int k = lo.z; k <= hi.z; ++k)
        {
            for     (int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    AllPrintToFile(pltfilefinal) << i << " " << Ubox(i,j,k,var,mat) << std::endl;
                }
            }
        }
    }
}

void PrintAllVarsTo1DGnuplotFile(CellArray& U, int picture, std::string filename)
{

    PrintTo1DGnuplotFile(U, filename,"/rho"     ,picture, RHO);
    PrintTo1DGnuplotFile(U, filename,"/p"       ,picture, P);
    PrintTo1DGnuplotFile(U, filename,"/u"       ,picture, VELOCITY);
    PrintTo1DGnuplotFile(U, filename,"/E"       ,picture, TOTAL_E);
    PrintTo1DGnuplotFile(U, filename,"/rhoU"    ,picture, RHOU);
    PrintTo1DGnuplotFile(U, filename,"/alpha0"  ,picture, ALPHA,0);

}

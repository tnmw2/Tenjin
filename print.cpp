#include "simulationheader.h"

void PrintTo1DGnuplotFile(CellArray& U, std::string filename, Variable var)
{
    for( MFIter mfi(U.data); mfi.isValid(); ++mfi )
    {
        const Box& box = mfi.validbox();

        FArrayBox& fab = U.data[mfi];

        Array4<Real> const& prop_arr = fab.array();

        const auto lo = lbound(box);
        const auto hi = ubound(box);

        BoxAccessCellArray Ubox(box,fab,prop_arr,U);


        for    		(int k = lo.z; k <= hi.z; ++k)
        {
            for     (int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    AllPrintToFile(filename) << i << " " << Ubox(i,j,k,var) << std::endl;
                }
            }
        }
    }
}

void PrintAllVarsTo1DGnuplotFile(CellArray& U, int picture, std::string filename)
{

    std::string pltfilevarstub = filename;
    std::string pltfilevar;

    pltfilevar = pltfilevarstub.append("/rho");
    const std::string pltfilerho = amrex::Concatenate(pltfilevar,picture,2);
    PrintTo1DGnuplotFile(U, pltfilerho, RHO);

    pltfilevarstub = filename;
    pltfilevar = pltfilevarstub.append("/u");
    const std::string pltfileu = amrex::Concatenate(pltfilevar,picture,2);
    PrintTo1DGnuplotFile(U, pltfileu, VELOCITY);

    pltfilevarstub = filename;
    pltfilevar = pltfilevarstub.append("/p");
    const std::string pltfilep = amrex::Concatenate(pltfilevar,picture,2);
    PrintTo1DGnuplotFile(U, pltfilep, P);

    pltfilevarstub = filename;
    pltfilevar = pltfilevarstub.append("/E");
    const std::string pltfileE = amrex::Concatenate(pltfilevar,picture,2);
    PrintTo1DGnuplotFile(U, pltfileE, TOTAL_E);

    pltfilevarstub = filename;
    pltfilevar = pltfilevarstub.append("/rhoU");
    const std::string pltfilerhoU = amrex::Concatenate(pltfilevar,picture,2);
    PrintTo1DGnuplotFile(U, pltfilerhoU, RHOU);

}

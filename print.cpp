#include "simulationheader.h"

void PrintTo1DGnuplotFile(MultiFab& phi, std::string filename, Variable var)
{
    for( MFIter mfi(phi); mfi.isValid(); ++mfi )
    {
        const Box& box = mfi.validbox();

        FArrayBox& fab = phi[mfi];

        Array4<Real> const& prop_arr = fab.array();

        const auto lo = lbound(box);
        const auto hi = ubound(box);


        for    		(int k = lo.z; k <= hi.z; ++k)
        {
            for     (int j = lo.y; j <= hi.y; ++j)
            {
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    AllPrintToFile(filename) << i << " " << prop_arr(i,j,k,var) << std::endl;
                }
            }
        }
    }
}

void PrintAllVarsTo1DGnuplotFile(MultiFab& phi, int picture, std::string filename)
{

    std::string pltfilevarstub = filename;
    std::string pltfilevar;

    pltfilevar = pltfilevarstub.append("/rho");
    const std::string pltfilerho = amrex::Concatenate(pltfilevar,picture,2);
    PrintTo1DGnuplotFile(phi, pltfilerho, RHO);

    pltfilevarstub = filename;
    pltfilevar = pltfilevarstub.append("/u");
    const std::string pltfileu = amrex::Concatenate(pltfilevar,picture,2);
    PrintTo1DGnuplotFile(phi, pltfileu, VELOCITY);

    pltfilevarstub = filename;
    pltfilevar = pltfilevarstub.append("/p");
    const std::string pltfilep = amrex::Concatenate(pltfilevar,picture,2);
    PrintTo1DGnuplotFile(phi, pltfilep, P);

}

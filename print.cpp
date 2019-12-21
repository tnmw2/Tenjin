#include "simulationheader.h"

void PrintTo1DGnuplotFile(CellArray& U, std::string filename, std::string name, int picture,  const Real* dx, const Real* prob_lo, Variable var, int mat =0, int row=0, int col=0)
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

        Real x,y,z;

        for 		(int k = lo.z; k <= hi.z; ++k)
        {
                    z = prob_lo[2] + (Real(k)+0.5)*dx[2];

            for 	(int j = lo.y; j <= hi.y; ++j)
            {
                    y = prob_lo[1] + (Real(j)+0.5)*dx[1];

                for (int i = lo.x; i <= hi.x; ++i)
                {
                    x = prob_lo[0] + (Real(i)+0.5)*dx[0];

                    AllPrintToFile(pltfilefinal) << x << " " << Ubox(i,j,k,var,mat,row,col) << std::endl;
                }
            }
        }
    }
}

void PrintAllVarsTo1DGnuplotFile(CellArray& U, int picture, std::string filename, const Real* dx, const Real* prob_lo)
{

    PrintTo1DGnuplotFile(U, filename,"/rho"     ,picture, dx, prob_lo, RHO);
    PrintTo1DGnuplotFile(U, filename,"/p"       ,picture, dx, prob_lo, P);
    PrintTo1DGnuplotFile(U, filename,"/u"       ,picture, dx, prob_lo, VELOCITY,0,0);
    PrintTo1DGnuplotFile(U, filename,"/v"       ,picture, dx, prob_lo, VELOCITY,0,1);
    PrintTo1DGnuplotFile(U, filename,"/w"       ,picture, dx, prob_lo, VELOCITY,0,2);
    PrintTo1DGnuplotFile(U, filename,"/sigmaxx" ,picture, dx, prob_lo, SIGMA,   0,0,0);
    PrintTo1DGnuplotFile(U, filename,"/sigmaxy" ,picture, dx, prob_lo, SIGMA,   0,0,1);
    PrintTo1DGnuplotFile(U, filename,"/sigmaxz" ,picture, dx, prob_lo, SIGMA,   0,0,2);

}

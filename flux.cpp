#include "simulationheader.h"

/** Takes the dot product of the velocity with a column of the stress tensor for a Cell.
 */
Real vdotsigma(Cell& U, int d, int m)
{
    Real result=0.0;

    for(int i=0;i<U.numberOfComponents;i++)
    {
        result+=U(VELOCITY,m,i)*U(SIGMA,m,i,d);
    }

    return result;
}

/** Takes the dot product of the velocity with a column of the stress tensor for a Box.
 */
Real vdotsigma(BoxAccessCellArray& U, int i, int j, int k, int d, int m)
{
    Real result=0.0;

    for(int row=0;row<U.numberOfComponents;row++)
    {
        result+=U(i,j,k,VELOCITY,m,row)*U(i,j,k,SIGMA,m,row,d);
    }

    return result;
}

/** Returns the equation flux of the conservative variables.
 */
Real flux(MaterialSpecifier n, Cell& U, Direction_enum d)
{
    switch(n.var)
    {
        case RHO:       	  return U(RHO,     n.mat)*            U(VELOCITY,n.mat,d);
        case RHOU: 			  return U(RHOU,    n.mat,n.row)*      U(VELOCITY,n.mat,d)-U(SIGMA,n.mat,n.row,d);
        case TOTAL_E:	   	  return U(TOTAL_E, n.mat)*            U(VELOCITY,n.mat,d)-vdotsigma(U,d,n.mat);
        case V_TENSOR:        return U(V_TENSOR,n.mat,n.row,n.col)*U(VELOCITY,n.mat,d)-U(V_TENSOR,n.mat,d,n.col)*U(VELOCITY,n.mat,n.row);
        default:   amrex::Print() << "Bad flux variable" << std::endl; exit(1);
    }

}

/** Calculates the Geometric flux for various different variables.
 */
Real geometricFlux(BoxAccessCellArray& U, int i, int j, int k, MaterialSpecifier n)
{
    switch(n.var)
    {
    case ALPHA:          return 0.0;
    case ALPHARHO:  	 return U(i,j,k,n)*U(i,j,k,VELOCITY,0,x);
    case ALPHARHOLAMBDA: return U(i,j,k,n)*U(i,j,k,VELOCITY,0,x);
    case ALPHARHOEPSILON:return U(i,j,k,n)*U(i,j,k,VELOCITY,0,x);

    case RHOU:
        if(n.row==0){       return U(i,j,k,n)*U(i,j,k,VELOCITY,0,x) - U(i,j,k,SIGMA,0,0,0) + U(i,j,k,SIGMA,0,2,2);}
        else if(n.row==1){  return U(i,j,k,n)*U(i,j,k,VELOCITY,0,x) - U(i,j,k,SIGMA,0,0,1);}
        else if(n.row==2){  return -2.0*U(i,j,k,SIGMA,0,0,2);}
        else {Print() << "Bad radial flux variable" << std::endl;}

    //case TOTAL_E:	    return U(i,j,k,n)*U(i,j,k,VELOCITY,0,x) - vdotsigma(U,i,j,k,x);

    case V_TENSOR:
        if(n.row==0){       return  (1.0/3.0)*U(i,j,k,n)*U(i,j,k,VELOCITY,0,x);}
        else if(n.row==1){  return  (1.0/3.0)*U(i,j,k,n)*U(i,j,k,VELOCITY,0,x);}
        else if(n.row==2){  return -(2.0/3.0)*U(i,j,k,n)*U(i,j,k,VELOCITY,0,x);}
        else {Print() << "Bad radial flux variable" << std::endl;}

    default:   amrex::Print() << "Bad radial flux variable" << std::endl; exit(1);
    }

}

#include "simulationheader.h"

class variableMatrix
{
    public:

    variableMatrix(Vector<MaterialSpecifier>& variableList);

    void generateMap();

    Real& operator() (MaterialSpecifier& a, MaterialSpecifier& b);

    Vector<Real> data;
    Vector<MaterialSpecifier>& list;
    std::map<MaterialSpecifier,int> map;
    int numberOfVariables;

};

bool MaterialSpecifier::operator < (const MaterialSpecifier& b) const
{
    return (var < b.var) || ((var == b.var) && ((mat < b.mat) || ((mat == b.mat) && ( (row < b.row) || ( (row == b.row) && (col < b.col))))));
}

bool MaterialSpecifier::operator == (const MaterialSpecifier& b) const
{
    return ((var == b.var) && (mat == b.mat) && (row == b.row) && (col == b.col));
}

variableMatrix::variableMatrix(Vector<MaterialSpecifier> &variableList) : list(variableList)
{
    numberOfVariables = variableList.size();

    data.resize(numberOfVariables*numberOfVariables);

    generateMap();
}

void variableMatrix::generateMap()
{
    int counter = 0;

    for(auto n: list)
    {
        map.insert(std::pair<MaterialSpecifier,int>(n,counter));

        counter++;
    }
}

Real& variableMatrix::operator ()(MaterialSpecifier& a, MaterialSpecifier& b)
{
    return data[map[a]*numberOfVariables+map[b]];
}

Real componentOfPrimitiveMatrix(Cell& Ustart, Cell& Uend, Real s, MaterialSpecifier& a, MaterialSpecifier& b, Direction_enum d)
{
    switch(a.var)
    {
    case ALPHA:
        if((b.var == a.var) && (a.mat == b.mat))
        {
            return (1.0-s)*Ustart(VELOCITY,0,d)+s*Uend(VELOCITY,0,d);
        }
        else
        {
            return 0.0;
        }
        break;
    case V_TENSOR:
        if(a.row == d)
        {
            if(a == b)
            {
                return (1.0-s)*Ustart(VELOCITY,0,d)+s*Uend(VELOCITY,0,d);
            }
            else if((b.var == VELOCITY) && (b.row == a.row))
            {
                return -2.0/3.0*((1.0-s)*Ustart(V_TENSOR,0,d,a.col)+s*Uend(V_TENSOR,0,d,a.col));
            }
            else
            {
                return 0.0;
            }
        }
        else
        {
            if(a == b)
            {
                return (1.0-s)*Ustart(VELOCITY,0,d)+s*Uend(VELOCITY,0,d);
            }
            else if((b.var == VELOCITY) && (b.row == a.row))
            {
                return -1.0*((1.0-s)*Ustart(V_TENSOR,0,d,a.col)+s*Uend(V_TENSOR,0,d,a.col));
            }
            else if((b.var == VELOCITY) && (b.row == d))
            {
                return 1.0/3.0*((1.0-s)*Ustart(a)+s*Uend(a));
            }
            else
            {
                return 0.0;
            }
        }
        break;
    default: Abort("Incorrect variable in path con flux"); return 0.0;
    }
}

Real calculateAMatrixIntegration(Cell& Ustart, Cell& Uend, Direction_enum d, MaterialSpecifier& row, MaterialSpecifier& col)
{
    Real s;

    const int orderOfQuadratureRule = 1;

    static const Real weight  [orderOfQuadratureRule] = {2.0};
    static const Real abscissa[orderOfQuadratureRule] = {0.0};

    //static const Real weight  [orderOfQuadratureRule] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    //static const Real abscissa[orderOfQuadratureRule] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};

    //static const Real weight  [orderOfQuadratureRule] = {(322.0-13.0*sqrt(70.0))/900.0,(322.0+13.0*sqrt(70.0))/900.0, 128.0/225.0, (322.0+13.0*sqrt(70.0))/900.0,(322.0-13.0*sqrt(70.0))/900.0};
    //static const Real abscissa[orderOfQuadratureRule] = {-(1.0/3.0)*sqrt(5.0+2.0*sqrt(10.0/7.0)), -(1.0/3.0)*sqrt(5.0-2.0*sqrt(10.0/7.0)), 0.0, (1.0/3.0)*sqrt(5.0-2.0*sqrt(10.0/7.0)), (1.0/3.0)*sqrt(5.0+2.0*sqrt(10.0/7.0))};

    Real sum = 0.0;

    for(int n = 0; n< orderOfQuadratureRule; n++)
    {
        s = 0.5*abscissa[n]+0.5;

        sum += 0.5*weight[n]*componentOfPrimitiveMatrix(Ustart,Uend,s,row,col,d);
    }

    return sum;
}

void calculatePathConservativeFluxes(BoxAccessCellArray& fluxboxL, BoxAccessCellArray& fluxboxR, BoxAccessCellArray& Ubox, BoxAccessCellArray& ULbox, BoxAccessCellArray& URbox, BoxAccessCellArray& UStarbox, ParameterStruct& parameters, Direction_enum d, const Real *dx)
{
    const auto lo = lbound(Ubox.box);
    const auto hi = ubound(Ubox.box);

    int extra[3] = {0,0,0};

    extra[d]=1;

    Material_type phase = (parameters.SOLID == 1 ? solid : fluid);

    for 		(int k = lo.z; k <= hi.z+extra[z]; ++k)
    {
        for 	(int j = lo.y; j <= hi.y+extra[y]; ++j)
        {
            for (int i = lo.x; i <= hi.x+extra[x]; ++i)
            {
                Cell UL(Ubox,i-extra[x],j-extra[y],k-extra[z],phase);
                Cell UR(Ubox,i,j,k,phase);
                Cell UStar(UStarbox,i,j,k,phase);

                for(auto row: Ubox.accessPattern.primitiveVariables)
                {
                    if(row.var == ALPHA || row.var == V_TENSOR)
                    {
                        fluxboxL(i,j,k,row) = 0.0;
                        fluxboxR(i,j,k,row) = 0.0;

                        for(auto col: Ubox.accessPattern.primitiveVariables)
                        {
                            if(col.var == ALPHA || col.var == V_TENSOR || col.var == VELOCITY)
                            {
                                fluxboxL(i,j,k,row) += calculateAMatrixIntegration(UL,UStar,d,row,col)*(UStar(col)-UL(col));
                                fluxboxR(i,j,k,row) += calculateAMatrixIntegration(UStar,UR,d,row,col)*(UR(col)-UStar(col));
                            }
                        }
                    }
                }
            }
        }
    }
}

void PCupdate(BoxAccessCellArray& fluxboxL, BoxAccessCellArray& fluxboxR, BoxAccessCellArray& Ubox, BoxAccessCellArray& U1box, ParameterStruct& parameters, Direction_enum d, Real dt, const Real* dx)
{
    const auto lo = lbound(Ubox.box);
    const auto hi = ubound(Ubox.box);

    for    			(auto n : Ubox.accessPattern.primitiveVariables)
    {
        if(n.var == ALPHA || n.var == V_TENSOR)
        {
            for 		(int k = lo.z; k <= hi.z; ++k)
            {
                for 	(int j = lo.y; j <= hi.y; ++j)
                {
                    for (int i = lo.x; i <= hi.x; ++i)
                    {
                        U1box(i,j,k,n) -= (dt/dx[d])*(fluxboxR(i,j,k,n) + fluxboxL.right(d,i,j,k,n));
                    }
                }
            }
        }
    }
}


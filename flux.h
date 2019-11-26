#ifndef FLUX
#define FLUX

#include "cell.h"
#include "cellarray.h"
#include "structdefinitions.h"


Real vdotsigma(Cell& U, int d, int m);
Real vdotsigma(BoxAccessCellArray& U, int i, int j, int k, int d, int m);
Real flux(MaterialSpecifier n, Cell& U, Direction_enum d);
Real geometricFlux(BoxAccessCellArray& U, int i, int j, int k, MaterialSpecifier n);



#endif // FLUX

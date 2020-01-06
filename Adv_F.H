#ifndef _Adv_F_H_
#define _Adv_F_H_

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif

  //void get_tagging_params(const int* name, const int* namlen);

  void phifill(BL_FORT_FAB_ARG_3D(state),
        const int* dlo, const int* dhi,
        const amrex_real* dx, const amrex_real* glo,
        const amrex_real* time, const int* bc);


#ifdef __cplusplus
}
#endif

#endif

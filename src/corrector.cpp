#include "flux.h"

#include <iostream>

void corrector_flux::solve(
        const state &le, const state &ri,
        const slope &ls, const slope &cs, const slope &rs,
        dir::Direction dir, const gasinfo &gas)
{
    (void)le;
    (void)ri;
    (void)ls;
    (void)cs;
    (void)rs;
    (void)dir;
    (void)gas;

    Vec F;
    F.fill(0);
    fden[0] = F[0];
    for (int i = 1; i < nc; i++) {
        fden[i] = F[i];
        fden[0] -= F[i];
    }
    fmom.x = F[nc + 0];
    fmom.y = F[nc + 1];
    fmom.z = F[nc + 2];
    fener = F[nc + 3];
}

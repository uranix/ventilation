#include "../include/flux.h"

#include <iostream>

struct EntropyFix {
    const double lcut;
    const double eps;
    EntropyFix(const double c, const double eps = .05) : lcut(c * eps), eps(eps) {  }
    double operator()(const double &lam) const {
        if (fabs(lam) < lcut)
            return (lcut + lcut + lam + lam) / (2 * lcut);
        return lam;
    }
};

void predictor_flux::solve(
        const state &le, const state &ri,
        dir::Direction dir, const gasinfo &gas)
{
    Vec UL, UR;
    Vec FL, FR, F;

    UL = slope::stateToVec(le);
    UR = slope::stateToVec(ri);

    FL = slope::stateToFlux(le, dir, gas);
    FR = slope::stateToFlux(ri, dir, gas);

    slope mid(le, ri, dir, gas);

    const auto &Om = mid.Omega();
    const auto &iOm = mid.iOmega();
#ifndef NDEBUG
    const auto &Id = Om * iOm;
    const auto &Id2 = Mat::Identity();
    if (!Id.isApprox(Id2, 1e-6)) {
        std::cout << "W = " << std::endl << Om << std::endl;
        std::cout << "iW = " << std::endl << iOm << std::endl;
        std::cout << "W . iW = " << std::endl << Id << std::endl;
        abort();
    }
#endif
    const auto &lam = mid.lambda();
    const double c = mid.c;
    F = .5 * (FL + FR + iOm * lam.cwiseAbs().unaryExpr(EntropyFix(c)).cwiseProduct(Om * (UL - UR)));

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

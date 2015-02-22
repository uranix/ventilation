#include "../include/flux.h"

#include <iostream>

struct LaxWendroff {
    const double dt_h;
    LaxWendroff(const double dt_h) : dt_h(dt_h) { }
    double operator()(const double lam) const {
        return 0.5 * (fabs(lam) - dt_h * lam * lam);
    }
};

struct LimiterKholodov4 {
    double chop(double v, double lo, double hi) const {
        if (lo > hi)
            std::swap(lo, hi);
        if (v > hi)
            return hi;
        if (v < lo)
            return lo;
        return v;
    }
    double alpha(double s) const {
        return (s + 1) * (s + 2) / 12;
    }
    double beta(double s) const {
        return (2 - s) * (s + 3) / 6;
    }
    double gamma(double s) const {
        return (s - 2) * (s + 1) / 12;
    }
    double slope(double dwl, double dwc, double dwr, double sl, double sc, double sr) {
        double v = dwl * alpha(sl) + dwc * beta(sc) + dwr * gamma(sr);
        double s = sc + 1e-6;
        v = chop(v, 0, 2 * dwl / s);
        v = chop(v, 0, 2 * dwc / (1 - s));
        return v;
    }
};

void corrector_flux::solve(
        const state &, const state &,
        const slope &ls, const slope &cs, const slope &rs,
        dir::Direction, const gasinfo &)
{
    LimiterKholodov4 L;

    const auto &mid = cs;
    const auto &iOm = mid.iOmega();
    Vec S;
    const auto &lam = mid.lambda();
    for (int i = 0; i < S.size(); i++)
        if (lam[i] > 0) {
            S[i] = L.slope(ls.dW[i], cs.dW[i], rs.dW[i],
                    ls.lambda(i) * dt_h,
                    lam[i] * dt_h,
                    rs.lambda(i) * dt_h);
        } else {
            S[i] = L.slope(rs.dW[i], cs.dW[i], ls.dW[i],
                    -rs.lambda(i) * dt_h,
                    -lam[i] * dt_h,
                    -ls.lambda(i) * dt_h);
        }
    const auto &R = lam.unaryExpr(LaxWendroff(dt_h)).cwiseProduct(S);
    const auto &F = iOm * R;

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

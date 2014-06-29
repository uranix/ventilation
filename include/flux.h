#ifndef __FLUX_H__
#define __FLUX_H__

#include "state.h"

struct avg_params {
    double density;
    vec velocity;
    double specific_energy;
    double pressure;

    double solve(const avg_params &left, const avg_params &right, const vec &norm);
};

template<int nc>
struct flux {
    double fdens[nc];
    vec fmom;
    double fener;

    double vmax;

    flux() {
        zero();
    }

    void zero() {
        for (int i = 0; i < nc; i++)
            fdens[i] = 0;
        fmom = vec(0);
        fener = 0;
        vmax = 0;
    }

    void add(const state<nc> &left, const state<nc> &right, const vec &norm, const double Sfrac, const gasinfo<nc> &gas) {
        avg_params la, ra, iface;

        la.density = left.density();
        la.velocity = left.velocity();
        la.specific_energy = left.specific_energy();
        la.pressure = gas.pressure(left);

        ra.density = right.density();
        ra.velocity = right.velocity();
        ra.specific_energy = right.specific_energy();
        ra.pressure = gas.pressure(right);

        double _vmax = iface.solve(la, ra, norm);
        if (_vmax > vmax)
            vmax = _vmax;

        double vn = iface.velocity.dot(norm);

        double theta[nc];

        if (vn > 0) {
            for (int i = 0; i < nc; i++)
                theta[i] = left.rho[i] / la.density;
        } else {
            for (int i = 0; i < nc; i++)
                theta[i] = right.rho[i] / ra.density;
        }

        for (int i = 0; i < nc; i++)
            fdens[i] += Sfrac * iface.density * theta[i] * vn;

        fmom += Sfrac * (iface.density * vn * iface.velocity + iface.pressure * norm);
        fener += Sfrac * vn * (
                iface.density * (iface.specific_energy + .5 * iface.velocity.norm2()) + iface.pressure
            );
    }

    void add_reflect(const state<nc> &inner, bool inner_is_left, const vec &norm, const double Sfrac, const gasinfo<nc> &gas) {
        state<nc> outer = inner;
        outer.rhou = inner.rhou - 2 * norm * dot(norm, inner.rhou);
        if (inner_is_left)
            add(inner, outer, norm, Sfrac, gas, true);
        else
            add(outer, inner, norm, Sfrac, gas, true);
    }
};

#endif

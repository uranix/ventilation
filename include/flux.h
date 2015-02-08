#ifndef __FLUX_H__
#define __FLUX_H__

#include "state.h"

/*
 * Gravity compensation reconstruction. Doesn't work, stub.
 * */
struct rec_params {
    double density;
    vec velocity;
    double specific_energy;
    double pressure;

    template<int nc>
    void reconstruct(const state<nc> &state, const gasinfo<nc> &gas, const vec &gh, const vec &norm) {
        (void)gh;
        (void)norm;
        density = state.density();
        velocity = state.velocity();
        specific_energy = state.specific_energy();
        pressure = gas.pressure(state);
    }
};

struct interface_flux {
    double fden;
    vec fmom;
    double fener;

    void solve(const rec_params &left, const rec_params &right, const vec &norm);
};

template<int nc>
struct flux {
    double fdens[nc];
    vec fmom;
    double fener;

    flux() {
        zero();
    }

    void zero() {
        for (int i = 0; i < nc; i++)
            fdens[i] = 0;
        fmom = vec(0);
        fener = 0;
    }

    void add_kernel(const double tl[], const double tr[], const rec_params &la, const rec_params &ra, const vec &norm, const double Sfrac) {
        interface_flux iface;

        iface.solve(la, ra, norm);

        const double rhovn = iface.fden;
        const double *theta = rhovn > 0 ? tl : tr;

        for (int i = 0; i < nc; i++)
            fdens[i] += Sfrac * iface.fden * theta[i];

        fmom += Sfrac * iface.fmom;
        fener += Sfrac * iface.fener;
    }

    void add(const state<nc> &left, const state<nc> &right, const vec &norm, const double Sfrac, const gasinfo<nc> &gas, const vec &gh) {
        rec_params la, ra;
        double tl[nc];
        double tr[nc];

        left.fractions(tl);
        right.fractions(tr);

        la.reconstruct(left, gas, gh, norm);
        ra.reconstruct(right, gas, -gh, norm);

        add_kernel(tl, tr, la, ra, norm, Sfrac);
    }

    void add_reflect(const state<nc> &inner, bool inner_is_left, const vec &norm, const double Sfrac, const gasinfo<nc> &gas, const vec &gh) {
        rec_params la, ra;
        double th[nc];

        inner.fractions(th);

        if (inner_is_left) {
            la.reconstruct(inner, gas, gh, norm);
            ra = la;
            ra.velocity = la.velocity - 2 * norm * dot(norm, la.velocity);
            add_kernel(th, th, la, ra, norm, Sfrac);
        } else {
            ra.reconstruct(inner, gas, -gh, norm);
            la = ra;
            la.velocity = ra.velocity - 2 * norm * dot(norm, ra.velocity);
            add_kernel(th, th, la, ra, norm, Sfrac);
        }
    }
};

#endif

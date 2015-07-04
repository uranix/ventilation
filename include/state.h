#ifndef __STATE_H__
#define __STATE_H__

#include "vec.h"
#include <vector>
#include <cassert>

constexpr int nc = NC;

/*
 * Cell state, stores conservative values.
 * */

struct gasinfo;

struct state {
    double rho[nc];
    vec rhou;
    double e;
#if TURBULENCE
    vec gradu[3]; // gradu[i] = du / dx_i
    double rhok, rhoeps;
    double Pk;
#endif

    state() {
        zero();
    }

    void zero() {
        rhou = vec(0);
        e = 0;
#if TURBULENCE
        rhok = 0;
        rhoeps = 0;
        Pk = 0;
#endif
        for (int i = 0; i < nc; i++)
            rho[i] = 0;
    }

    void from_rue(const std::vector<double> &r, const vec &u, double eps);
    void from_rup(const std::vector<double> &r, const vec &u, double p, const gasinfo &gas);
    void from_ruT(const std::vector<double> &r, const vec &u, double T, const gasinfo &gas);

    double density() const {
        double sum = 0;
        for (int i = 0; i < nc; i++)
            sum += rho[i];
        return sum;
    }

    void fractions(double theta[]) const {
        double dens = density();
        for (int i = 0; i < nc; i++)
            theta[i] = rho[i] / dens;
    }

    vec velocity() const {
        double r = density();
        return rhou / r;
    }

    double specific_energy() const {
        double r = density();
        return (e - 0.5 * rhou.norm2() / r) / r;
    }
#if TURBULENCE
    double turb_energy() const {
        return rhok / density();
    }
    double turb_dissipation() const {
        return rhoeps / density();
    }
    double turb_viscosity() const;
#endif
};

template<class T>
class optional {
    T *ptr;
public:
    optional(std::nullptr_t) {
        ptr = nullptr;
    }
    optional(T &ref) {
        ptr = &ref;
    }
    operator bool() const {
        return ptr != nullptr;
    }
    bool operator!() const {
        return ptr == nullptr;
    }
    T &operator*() {
        assert(!!*this);
        return *ptr;
    }
};

#endif

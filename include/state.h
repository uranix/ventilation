#ifndef __STATE_H__
#define __STATE_H__

#include "vec.h"
#include <vector>

/*
 * Cell state, stores conservative values.
 * */

template<int nc>
struct gasinfo;

template<int nc>
struct state {
    double rho[nc];
    vec rhou;
    double rhoE;

    state() {
        zero();
    }

    void zero() {
        rhou = vec(0);
        rhoE = 0;
        for (int i = 0; i < nc; i++)
            rho[i] = 0;
    }

    void from_rue(const std::vector<double> &r, const vec &u, double eps) {
        double rs = 0;
        for (int i = 0; i < nc; i++) {
            rho[i] = r[i];
            rs += rho[i];
        }

        rhou = rs * u;
        rhoE = rs * (eps + 0.5 * u.norm2());
    }

    void from_ruT(const std::vector<double> &r, const vec &u, double T, const gasinfo<nc> &gas);

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
        return (rhoE - 0.5 * rhou.norm2() / r) / r;
    }
};

template<int nc>
struct gasinfo {
    double Rspecific[nc]; /* J / kg / K */
    double gamma[nc];
    double visc[nc];

    gasinfo() {
    }

    void set(int i, double molar_mass, double gamma_factor, double viscosity) {
        const double Runi = 8314.4;

        Rspecific[i] = Runi / molar_mass;
        gamma[i] = gamma_factor;
        visc[i] = viscosity;
    }

    double density(const state<nc> &st) const {
        return st.density();
    }

    double specific_energy(const state<nc> &st) const {
        return st.specific_energy();
    }

    double pressure(const state<nc> &st) const {
        return (gamma_factor(st) - 1) * st.density() * st.specific_energy();
    }

    double gamma_factor(const state<nc> &st) const {
        double Cv = 0, Cp = 0;

        for (int i = 0; i < nc; i++) {
            double x = st.rho[i] * Rspecific[i] / (gamma[i] - 1);
            Cv += x;
            Cp += x * gamma[i];
        }

        return Cp / Cv;
    }

    double heat_capacity_volume(const state<nc> &st) const {
        double Cv = 0;

        for (int i = 0; i < nc; i++) {
            double x = st.rho[i] * Rspecific[i] / (gamma[i] - 1);
            Cv += x;
        }

        return Cv / st.density();
    }

    double viscosity(const state<nc> &st) const {
        /* Herning-Zipprer */
        double x[nc], xsum = 0;
        for (int i = 0; i < nc; i++) {
            x[i] = st.rho[i] * Rspecific[i];
            xsum += x[i];
        }
        for (int i = 0; i < nc; i++) {
            x[i] /= xsum;
        }
        double mumix = 0;
        for (int i = 0; i < nc; i++) {
            double denom = 0;
            for (int j = 0; j < nc; j++)
                denom += x[j] * sqrt(visc[i] / visc[j]);
            mumix += x[i] * visc[i] / denom;
        }

        return mumix;
    }

    double temperature(const state<nc> &st) const {
        return st.specific_energy() / heat_capacity_volume(st);
    }

    double sound_speed(const state<nc> &st) const {
        const double g = gamma_factor(st);
        return sqrt(g * (g - 1) * specific_energy(st));
    }
};

template<int nc>
void state<nc>::from_ruT(const std::vector<double> &r, const vec &u, double T, const gasinfo<nc> &gas) {
    from_rue(r, u, 0);

    double eps = gas.heat_capacity_volume(*this) * T;

    from_rue(r, u, eps);
}

#endif

#ifndef __GASINFO_H__
#define __GASINFO_H__

#include "../include/state.h"

struct gasinfo {
    double Rspecific[nc]; /* J / kg / K */
    double gamma[nc];
    double molar[nc];
    double beta[nc];
    double visc[nc];

    gasinfo() {
    }

    void set(int i, double molar_mass, double gamma_ratio, double viscosity);
    double density(const state &st) const;
    double specific_energy(const state &st) const;
    double pressure(const state &st) const;
    double gamma_ratio(const state &st) const;
    double heat_capacity_volume(const state &st) const;
    double viscosity(const state &st) const;
    double temperature(const state &st) const;
    double beta_ratio(const state &st) const;
    double dbetadtheta(const state &st, int i) const;
    double molar_mass(const state &st) const;
    double sound_speed(const state &st) const;
};

#endif

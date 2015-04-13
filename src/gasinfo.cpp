#include "../include/gasinfo.h"

void gasinfo::set(int i, double molar_mass, double gamma_ratio, double viscosity) {
    const double Runi = 8314.4;

    Rspecific[i] = Runi / molar_mass;
    gamma[i] = gamma_ratio;
    molar[i] = molar_mass;
    beta[i] = 1 / (gamma_ratio - 1);
    visc[i] = viscosity;
}

double gasinfo::density(const state &st) const {
    return st.density();
}

double gasinfo::specific_energy(const state &st) const {
    return st.specific_energy();
}

double gasinfo::pressure(const state &st) const {
    return (gamma_ratio(st) - 1) * st.density() * st.specific_energy();
}

double gasinfo::gamma_ratio(const state &st) const {
    double Cv = 0, Cp = 0;

    for (int i = 0; i < nc; i++) {
        double x = st.rho[i] * Rspecific[i];
        Cv += x * beta[i];
        Cp += x * (beta[i] + 1);
    }

    return Cp / Cv;
}

double gasinfo::heat_capacity_volume(const state &st) const {
    double Cv = 0;

    for (int i = 0; i < nc; i++) {
        double x = st.rho[i] * Rspecific[i] * beta[i];
        Cv += x;
    }

    return Cv / st.density();
}

double gasinfo::viscosity(const state &st) const {
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

double gasinfo::temperature(const state &st) const {
    return st.specific_energy() / heat_capacity_volume(st);
}

double gasinfo::beta_ratio(const state &st) const {
    double Cv = 0, MR = 0;

    for (int i = 0; i < nc; i++) {
        double x = st.rho[i] * Rspecific[i];
        Cv += x * beta[i];
        MR += x;
    }

    return Cv / MR;
}

double gasinfo::dbetadtheta(const state &st, int i) const {
    const double M = molar_mass(st);
    const double b = beta_ratio(st);
    return M * ((b - beta[0]) / molar[0] - (b - beta[i]) / molar[i]);
}

double gasinfo::molar_mass(const state &st) const {
    double rho_M = 0, rho = 0;

    for (int i = 0; i < nc; i++) {
        rho_M += st.rho[i] / molar[i];
        rho += st.rho[i];
    }

    return rho / rho_M;
}

double gasinfo::sound_speed(const state &st) const {
    const double g = gamma_ratio(st);
    return sqrt(g * (g - 1) * specific_energy(st));
}

#include "../include/state.h"
#include "../include/gasinfo.h"

void state::from_ruT(const std::vector<double> &r, const vec &u, double T, const gasinfo &gas) {
    from_rue(r, u, 0);

    double eps = gas.heat_capacity_volume(*this) * T;

    from_rue(r, u, eps);
}

void state::from_rup(const std::vector<double> &r, const vec &u, double p, const gasinfo &gas) {
    from_rue(r, u, 0);

    double eps = gas.beta_ratio(*this) * p / density();

    from_rue(r, u, eps);
}

#if TURBULENCE
double state::turb_viscosity() const {
    const double Cmu = 0.09;

    // For ghost cells
    if (rhok == 0 && rhoeps == 0)
        return -1;

    return Cmu * rhok * rhok / rhoeps;
}
#endif

void state::from_rue(const std::vector<double> &r, const vec &u, double eps) {
    double rs = 0;
    for (int i = 0; i < nc; i++) {
        rho[i] = r[i];
        rs += rho[i];
    }

    rhou = rs * u;
    e = rs * (eps + 0.5 * u.norm2());

#if TURBULENCE
    // Estimated for
    // L = 1m, v = 10 m/s, Re = 7.5e5
    const double U = 10;
    const double L = 1;
    const double nu = 1.75e-5;
    const double Re = density() * U / nu;
    const double I = 0.16 * std::pow(Re, -.125);
    const double Cmu = 0.09;
    const double ell = 0.07 * L;

    rhok = 1.5 * std::pow(U * I, 2) * rs;
    rhoeps = std::pow(Cmu, .75) * std::pow(rhok / rs, 1.5) / ell * rs;
#endif
}

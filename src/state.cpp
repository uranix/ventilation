#include "../include/state.h"
#include "../include/gasinfo.h"

void state::from_ruT(const std::vector<double> &r, const vec &u, double T, const gasinfo &gas) {
    from_rue(r, u, 0);

    double eps = gas.heat_capacity_volume(*this) * T;

    from_rue(r, u, eps);
}

void state::from_rue(const std::vector<double> &r, const vec &u, double eps) {
    double rs = 0;
    for (int i = 0; i < nc; i++) {
        rho[i] = r[i];
        rs += rho[i];
    }

    rhou = rs * u;
    e = rs * (eps + 0.5 * u.norm2());
}

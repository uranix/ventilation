#include "flux.h"

#include "riemann/euler.h"

using namespace riemann_solver;

double avg_params::solve(const avg_params &left, const avg_params &right, const vec &norm) {
    typedef euler<double, 3>::vec rvec;
    rvec v1, v2, n;

    v1[0] = left.velocity.x;
    v1[1] = left.velocity.y;
    v1[2] = left.velocity.z;

    v2[0] = right.velocity.x;
    v2[1] = right.velocity.y;
    v2[2] = right.velocity.z;

    n[0] = norm.x;
    n[1] = norm.y;
    n[2] = norm.z;

    euler<double, 3> solver(
        left.density, v1, left.pressure, 1 + left.pressure / left.density / left.specific_energy,
        right.density, v2, right.pressure, 1 + right.pressure / right.density / right.specific_energy,
        n);

    double r, p, eps, xi = 0;
    solver.get(1, &xi, &r, &v1, &eps, &p);

    density = r;
    pressure = p;
    specific_energy = eps;

    velocity.x = v1[0];
    velocity.y = v1[1];
    velocity.z = v1[2];

    return solver.get_1d_solver().max_speed();
}

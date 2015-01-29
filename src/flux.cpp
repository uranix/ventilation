#include "flux.h"

#include "riemann/euler.h"
#include <iostream>

using namespace riemann_solver;

#if PRECISE_RIEMANN
double avg_flux::solve(const avg_params &left, const avg_params &right, const vec &norm) {
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

    double tol = 1e-6;
    const int maxit = 5;
    euler<double, 3> solver(
        left.density, v1, left.pressure, 1 + left.pressure / left.density / left.specific_energy,
        right.density, v2, right.pressure, 1 + right.pressure / right.density / right.specific_energy,
        n, tol, maxit);

    double r, p, eps, xi = 0;
    solver.get(1, &xi, &r, &v1, &eps, &p);

    double vn = v1[0] * n[0] + v1[1] * n[1] + v1[2] * n[2];
    double vv = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];

    fden = r * vn;

    fmom.x = fden * v1[0] + p * n[0];
    fmom.y = fden * v1[1] + p * n[1];
    fmom.z = fden * v1[2] + p * n[2];

    fener = fden * (eps + 0.5 * vv) + p * vn;

    return solver.get_1d_solver().max_speed();
}
#else
double avg_flux::solve(const avg_params &left, const avg_params &right, const vec &norm) {
    double gl = left.pressure / left.density / left.specific_energy + 1;
    double gr = right.pressure / right.density / right.specific_energy + 1;
    double cleft = sqrt(gl * left.pressure / left.density);
    double cright = sqrt(gr * right.pressure / right.density);

    double vnl = left.velocity.dot(norm);
    double vnr = right.velocity.dot(norm);
    double aleft = cleft + fabs(vnl);
    double aright = cright + fabs(vnr);

    double amax = std::max(aleft, aright);

    double v2l = left.velocity.dot(left.velocity);
    double v2r = right.velocity.dot(left.velocity);
/*
    double fdenl = left.density * vnl;
    vec fmoml = left.density * left.velocity * vnl + left.pressure * norm;
    double fenerl = (left.density * (0.5 * v2l + left.specific_energy) + left.pressure) * vnl;

    double fdenr = right.density * vnr;
    vec fmomr = right.density * right.velocity * vnr + right.pressure * norm;
    double fenerr = (right.density * (0.5 * v2r + right.specific_energy) + right.pressure) * vnr;
*/
    fden = 0.5 * (left.density * vnl + right.density * vnr
                + amax * (left.density - right.density));
    fmom = 0.5 * (left.density * left.velocity * vnl + right.density * right.velocity * vnr
                + left.pressure * norm + right.pressure * norm
                + amax * (left.density * left.velocity - right.density * right.velocity));
    fener = 0.5 * (left.density * (left.specific_energy + 0.5 * v2l) * vnl + right.density * (right.specific_energy + 0.5 * v2r) * vnr
                + left.pressure * vnl + right.pressure * vnr
                + amax * (left.density * (left.specific_energy + 0.5 * v2l) - right.density * (right.specific_energy + 0.5 * v2r)));
    return amax;
}
#endif

#include "riemann/euler.h"

using namespace riemann_solver;

template<int nc>
void solver_flux<nc>::solve(const state<nc> &le, const state<nc> &ri, const vec &norm, const gasinfo<nc> &gas) {
    typedef euler<double, 3>::vec rvec;
    rvec v1, v2, n;

    v1[0] = le.velocity().x;
    v1[1] = le.velocity().y;
    v1[2] = le.velocity().z;

    v2[0] = ri.velocity().x;
    v2[1] = ri.velocity().y;
    v2[2] = ri.velocity().z;

    n[0] = norm.x;
    n[1] = norm.y;
    n[2] = norm.z;

    double tol = 1e-6;
    const int maxit = 5;
    euler<double, 3> solver(
        le.density(), v1, gas.pressure(le), gas.gamma_factor(le),
        ri.density(), v2, gas.pressure(ri), gas.gamma_factor(ri),
        n, tol, maxit);

    double r, p, eps, xi = 0;
    solver.get(1, &xi, &r, &v1, &eps, &p);

    double vn = v1[0] * n[0] + v1[1] * n[1] + v1[2] * n[2];
    double vv = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];

    double fdentot = r * vn;

    fmom.x = fdentot * v1[0] + p * n[0];
    fmom.y = fdentot * v1[1] + p * n[1];
    fmom.z = fdentot * v1[2] + p * n[2];

    fener = fdentot * (eps + 0.5 * vv) + p * vn;

    double t[nc];
    if (vn > 0)
        le.fractions(t);
    else
        ri.fractions(t);
    for (int i = 0; i < nc; i++)
        fden[i] = fdentot * t[i];
}

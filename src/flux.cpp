#include "flux.h"

#include "riemann/euler.h"
#include <Eigen/Core>
#include <iostream>

using namespace riemann_solver;

#if PRECISE_RIEMANN
void interface_flux::solve(const rec_params &le, const rec_params &ri, const vec &norm) {
    typedef euler<double, 3>::vec rvec;
    rvec v1, v2, n;

    v1[0] = le.velocity.x;
    v1[1] = le.velocity.y;
    v1[2] = le.velocity.z;

    v2[0] = ri.velocity.x;
    v2[1] = ri.velocity.y;
    v2[2] = ri.velocity.z;

    n[0] = norm.x;
    n[1] = norm.y;
    n[2] = norm.z;

    double tol = 1e-6;
    const int maxit = 5;
    euler<double, 3> solver(
        le.density, v1, le.pressure, 1 + le.pressure / le.density / le.specific_energy,
        ri.density, v2, ri.pressure, 1 + ri.pressure / ri.density / ri.specific_energy,
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
}
#else
void interface_flux::solve(const rec_params &le, const rec_params &ri, const vec &n) {
    Eigen::Matrix<double, 5, 1> UL, UR, U, lam;
    Eigen::Matrix<double, 5, 1> FL, FR, F;
    Eigen::Matrix<double, 5, 5> Om, iOm;

    double q2l = le.velocity.dot(le.velocity);
    double q2r = ri.velocity.dot(ri.velocity);
    double vnl = le.velocity.dot(n);
    double vnr = ri.velocity.dot(n);

    UL[0] = le.density;
    UL[1] = le.density * le.velocity.x;
    UL[2] = le.density * le.velocity.y;
    UL[3] = le.density * le.velocity.z;
    UL[4] = le.density * (le.specific_energy + .5 * q2l);

    UR[0] = ri.density;
    UR[1] = ri.density * ri.velocity.x;
    UR[2] = ri.density * ri.velocity.y;
    UR[3] = ri.density * ri.velocity.z;
    UR[4] = ri.density * (ri.specific_energy + .5 * q2r);

    FL[0] = UL[0] * vnl;
    FL[1] = UL[1] * vnl + le.pressure * n.x;
    FL[2] = UL[2] * vnl + le.pressure * n.y;
    FL[3] = UL[3] * vnl + le.pressure * n.z;
    FL[4] = UL[4] * vnl + le.pressure * vnl;

    FR[0] = UR[0] * vnr;
    FR[1] = UR[1] * vnr + ri.pressure * n.x;
    FR[2] = UR[2] * vnr + ri.pressure * n.y;
    FR[3] = UR[3] * vnr + ri.pressure * n.z;
    FR[4] = UR[4] * vnr + ri.pressure * vnr;

    U = 0.5 * (UL + UR);
    double rho = U[0];
    vec u(U[1], U[2], U[3]);
    u *= (1 / rho);
    double q2 = u.dot(u);
    double vn = u.dot(n);
    double eps = U[4];
    eps /= rho;
    eps -= .5 * q2;
    double kl = le.pressure / le.specific_energy / le.density;
    double kr = ri.pressure / ri.specific_energy / ri.density;
    //double k = 2 * kr * kl / (kr + kl);
    double k;
    if (vn > 0)
        k = kl;
    else
        k = kr;
    double c = sqrt(k * (k + 1) * eps);

    Om <<
        -u.x,                   1,                  0,                  0,                  0,
        -u.y,                   0,                  1,                  0,                  0,
        -u.z,                   0,                  0,                  1,                  0,
        .5 * q2 + c * vn / k,   -u.x - c/k * n.x,   -u.y - c/k * n.y,   -u.z - c/k * n.z,   1,
        .5 * q2 - c * vn / k,   -u.x + c/k * n.x,   -u.y + c/k * n.y,   -u.z + c/k * n.z,   1;

    int idx;
    if (n.x > .5)
        idx = 0;
    if (n.y > .5)
        idx = 1;
    if (n.z > .5)
        idx = 2;

    double s = c*c/k;

    Om.row(idx) << -s - q2/2 + vn*vn, -vn * n.x, -vn * n.y, -vn * n.z, 1;

    lam << vn, vn, vn, vn - c, vn + c;
//    const double amax = fabs(vn) + c;
//    lam << amax, amax, amax, amax, amax;

    iOm <<
        u.x,                u.y,                u.z,                .5,                     .5,
        u.x * u.x + s,      u.x * u.y,          u.x * u.z,          (u.x - c * n.x) / 2,    (u.x + c * n.x) / 2,
        u.x * u.y,          u.y * u.y + s,      u.y * u.z,          (u.y - c * n.y) / 2,    (u.y + c * n.y) / 2,
        u.x * u.z,          u.y * u.z,          u.z * u.z + s,      (u.z - c * n.z) / 2,    (u.z + c * n.z) / 2,
        (s + q2/2) * u.x,   (s + q2/2) * u.y,   (s + q2/2) * u.z,   q2/4 - c*vn/2 + s/2,    q2/4 + c*vn/2 + s/2;

    iOm.col(idx) << -1, -u.x, -u.y, -u.z, -.5*q2;

    iOm *= k / (c * c);

    F = .5 * (FL + FR + iOm * lam.cwiseAbs().cwiseProduct(Om * (UL - UR)));

    fden = F[0];
    fmom.x = F[1];
    fmom.y = F[2];
    fmom.z = F[3];
    fener = F[4];
}
#endif

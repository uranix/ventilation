#include "../include/fan.h"

namespace objects {

struct fan_solver {
    double P1, P2, U, R1, R2;

    double leftP(double rL, double uL, double pL, double U, double gL, double &R) {
        double Pleft;
        if (uL >= U) {
            double W = uL - U;
            Pleft = pL + 0.25 * (gL + 1) * rL * W * W + 0.25 * W *
                std::sqrt(16 * rL * pL * gL + std::pow((gL + 1) * W * rL, 2));
            double m = std::sqrt(0.5 * rL * ((gL + 1) * Pleft + (gL - 1) * pL));
            R = rL * m / (m - W * rL);
        } else {
            double W = U - uL;
            double c = sqrt(gL * pL / rL);
            Pleft = pL * std::pow(1 - 0.5 * (gL - 1) * W / c, 2 * gL / (gL - 1));
            R = rL * std::pow(Pleft / pL, 1 / gL);
        }
        return Pleft;
    }

    double rightP(double rR, double uR, double pR, double U, double gR, double &R) {
        double Pright;
        if (U >= uR) {
            double W = U - uR;
            Pright = pR + 0.25 * (gR + 1) * rR * W * W + 0.25 * W *
                std::sqrt(16 * rR * pR * gR + std::pow((gR + 1) * W * rR, 2));
            double m = std::sqrt(0.5 * rR * ((gR + 1) * Pright + (gR - 1) * pR));
            R = rR * m / (m - W * rR);
        } else {
            double W = uR - U;
            double c = sqrt(gR * pR / rR);
            Pright = pR * std::pow(1 - 0.5 * (gR - 1) * W / c, 2 * gR / (gR - 1));
            R = rR * std::pow(Pright / pR, 1 / gR);
        }
        return Pright;
    }

    fan_solver(
        double rl, double ul, double pl, double gl,
        double rr, double ur, double pr, double gr,
        double U)
    {
        this->U = U;
        P1 = leftP(rl, ul, pl, U, gl, R1);
        P2 = rightP(rr, ur, pr, U, gr, R2);
    }
};

void fan::compute_special_flux(dir::Direction dir, const double) {
    if (dir != this->dir)
        return;
    const auto &left = val(dir, 0);
    const auto &right = val(dir, 1);

    const auto &vl = left.velocity();
    const auto &vr = right.velocity();

    const auto &gas = this->gas();

    fan_solver fs(
            left.density(), vl(dir), gas.pressure(left), gas.gamma_ratio(left),
            right.density(), vr(dir), gas.pressure(right), gas.gamma_ratio(right),
            U);

    vec n(dir);
    double t[nc];
    double R, g;
    vec V;
    if (U > 0) {
        left.fractions(t);
        R = fs.R1;
        V = left.velocity();
        V(dir) = U;
        g = gas.gamma_ratio(left);
    } else {
        right.fractions(t);
        R = fs.R2;
        V = right.velocity();
        V(dir) = U;
        g = gas.gamma_ratio(right);
    }
    for (int i = 0; i < nc; i++)
        FL.fdens[i] = FR.fdens[i] = R * U * t[i];
    FL.fmom = R * U * V + fs.P1 * n;
    FR.fmom = R * U * V + fs.P2 * n;

    double cp = g / (g - 1);
    FL.fener = (cp * fs.P1 + 0.5 * R * V.norm2()) * U;
    FR.fener = (cp * fs.P2 + 0.5 * R * V.norm2()) * U;
}

void fan::integrate(state &cell, const flux &left, const flux &right, dir::Direction dir, double h, const double, const double dt) {
    if (&cell == &this->val(dir, 0)) {
        for (int i = 0; i < nc; i++)
            cell.rho[i] -= dt * (FL.fdens[i] - left.fdens[i]) / h;

        cell.rhou(dir) -= dt * (FL.fmom(dir) - left.fmom(dir)) / h;
        cell.e -= dt * (FL.fener - left.fener) / h;
        return;
    }
    if (&cell == &this->val(dir, 1)) {
        for (int i = 0; i < nc; i++)
            cell.rho[i] -= dt * (right.fdens[i] - FR.fdens[i]) / h;

        cell.rhou(dir) -= dt * (right.fmom(dir) - FR.fmom(dir)) / h;
        cell.e -= dt * (right.fener - FR.fener) / h;
        return;
    }
    throw;
}

}

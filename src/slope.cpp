#include "../include/slope.h"
#include "../include/gasinfo.h"

Vec slope::stateToVec(const state &st) {
    Vec U;
    U[0] = st.density();
    for (int i = 1; i < nc; i++)
        U[i] = st.rho[i];
    U[nc + 0] = st.rhou.x;
    U[nc + 1] = st.rhou.y;
    U[nc + 2] = st.rhou.z;
    U[nc + 3] = st.e;
    return U;
}

state slope::vecToState(const Vec &U) {
    state st;
    st.rho[0] = U[0];
    for (int i = 1; i < nc; i++) {
        st.rho[i] = U[i];
        st.rho[0] -= U[i];
    }
    st.rhou = vec(U[nc], U[nc + 1], U[nc + 2]);
    st.e = U[nc + 3];
    return st;
}

slope::slope() : c(1), v(0), dir(dir::X) {
    dW.fill(0);
}

slope::slope(
        const state &le, const state &ri,
        const dir::Direction dir, const gasinfo &gas)
    : dir(dir)
{
    const auto &UL = stateToVec(le);
    const auto &UR = stateToVec(ri);

    const auto &st = vecToState(0.5 * (UL + UR));
    st.fractions(theta);
    gamma = gas.gamma_factor(st);
    c = gas.sound_speed(st);
    v = st.velocity();

    for (int i = 1; i < nc; i++)
        dbdth[i] = gas.dbetadtheta(st, i);

    dW = Omega() * (UR - UL);
}

Mat slope::Omega() const {
    Mat Om;
    Om.fill(0);

    const vec n(dir);
    const double vn = v.dot(n);
    const double q2 = v.dot(v);
    const double g = gamma;
    const double b = 1 / (g - 1);
    const double c2b = c*c*b;

    for (int i = 1; i < nc; i++) {
        Om(i - 1, 0) = -theta[i];
        Om(i - 1, i) = 1;
    }

    Om(nc - 1, 0) = n.x > 0.5 ? -vn * vn + q2/2 + c2b : -v.x;
    Om(nc    , 0) = n.y > 0.5 ? -vn * vn + q2/2 + c2b : -v.y;
    Om(nc + 1, 0) = n.z > 0.5 ? -vn * vn + q2/2 + c2b : -v.z;
    Om(nc - 1, nc    ) = n.x > 0.5 ? vn : 1;
    Om(nc    , nc + 1) = n.y > 0.5 ? vn : 1;
    Om(nc + 1, nc + 2) = n.z > 0.5 ? vn : 1;
    Om(nc - 1, nc + 3) = -n.x;
    Om(nc    , nc + 3) = -n.y;
    Om(nc + 1, nc + 3) = -n.z;

    Om(nc + 2, 0) = q2 / 2 - c * vn * b;
    Om(nc + 3, 0) = q2 / 2 + c * vn * b;

    for (int i = 1; i < nc; i++) {
        Om(nc + 2, i) = -c * c / g * dbdth[i];
        Om(nc + 3, i) = -c * c / g * dbdth[i];
        Om(nc + 2, 0) += c * c / g * dbdth[i] * theta[i];
        Om(nc + 3, 0) += c * c / g * dbdth[i] * theta[i];
    }

    Om(nc + 2, nc + 0) = -v.x + n.x * c * b;
    Om(nc + 2, nc + 1) = -v.y + n.y * c * b;
    Om(nc + 2, nc + 2) = -v.z + n.z * c * b;

    Om(nc + 3, nc + 0) = -v.x - n.x * c * b;
    Om(nc + 3, nc + 1) = -v.y - n.y * c * b;
    Om(nc + 3, nc + 2) = -v.z - n.z * c * b;

    Om(nc + 2, nc + 3) = Om(nc + 3, nc + 3) = 1;

    return Om;
}

Mat slope::iOmega() const {
    Mat iOm;
    iOm.fill(0);

    const vec n(dir);
    const double vn = v.dot(n);
    const double q2 = v.dot(v);
    const double g = gamma;
    const double b = 1 / (g - 1);
    const double c2b = c*c*b;

    Vec col, row;

    col[0] = 1;
    for (int i = 1; i < nc; i++)
        col[i] = theta[i];
    col[nc + 0] = v.x;
    col[nc + 1] = v.y;
    col[nc + 2] = v.z;
    col[nc + 3] = q2/2 + c2b;

    for (int i = 1; i < nc; i++)
        row[i - 1] = c * c / g * dbdth[i];

    row[nc - 1] = n.x > .5 ? 1 : v.x;
    row[nc + 0] = n.y > .5 ? 1 : v.y;
    row[nc + 1] = n.z > .5 ? 1 : v.z;
    row[nc + 2] = .5;
    row[nc + 3] = .5;

    iOm = (col * row.transpose()) / c2b;
    for (int i = 1; i < nc; i++)
        iOm(i, i - 1) += 1;
    iOm(nc    , nc - 1) += 1 - n.x;
    iOm(nc + 1, nc    ) += 1 - n.y;
    iOm(nc + 2, nc + 1) += 1 - n.z;
    iOm(nc + 3, nc - 1) -= n.x;
    iOm(nc + 3, nc    ) -= n.y;
    iOm(nc + 3, nc + 1) -= n.z;

    const double z = .5 / c / b;
    iOm(nc    , nc + 2) += n.x * z;
    iOm(nc + 1, nc + 2) += n.y * z;
    iOm(nc + 2, nc + 2) += n.z * z;
    iOm(nc + 3, nc + 2) += vn * z;
    iOm(nc    , nc + 3) -= n.x * z;
    iOm(nc + 1, nc + 3) -= n.y * z;
    iOm(nc + 2, nc + 3) -= n.z * z;
    iOm(nc + 3, nc + 3) -= vn * z;

    return iOm;
}

Vec slope::lambda() const {
    const double vn = v(dir);

    Vec lam;
    for (int i = 0; i < nc + 2; i++)
        lam[i] = vn;

    lam[nc + 2] = vn + c;
    lam[nc + 3] = vn - c;

    return lam;
}

double slope::lambda(size_t i) const {
    const double vn = v(dir);
    if (i < nc + 2)
        return vn;
    if (i == nc + 2)
        return vn + c;
    return vn - c;
}

Vec slope::stateToFlux(const state &st, dir::Direction dir, const gasinfo &gas) {
    const vec n(dir);
    Vec F;
    double rho = st.density();
    double vn = st.rhou(dir) / rho;
    double p = gas.pressure(st);
    F[0] = st.rhou(dir);
    for (int i = 1; i < nc; i++)
        F[i] = F[0] * st.rho[i] / rho;
    F[nc + 0] = st.rhou.x * vn + p * n.x;
    F[nc + 1] = st.rhou.y * vn + p * n.y;
    F[nc + 2] = st.rhou.z * vn + p * n.z;
    F[nc + 3] = (st.e + p) * vn;
    return F;
}

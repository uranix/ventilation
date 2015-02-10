#include <Eigen/Core>
#include <Eigen/LU>

template<int nc> using Vec = Eigen::Matrix<double, nc + 4, 1>;
template<int nc> using Mat = Eigen::Matrix<double, nc + 4, nc + 4>;

template<int nc>
Vec<nc> formVector(const state<nc> &st) {
    Vec<nc> U;
    U[0] = st.density();
    for (int i = 1; i < nc; i++)
        U[i] = st.rho[i];
    U[nc + 0] = st.rhou.x;
    U[nc + 1] = st.rhou.y;
    U[nc + 2] = st.rhou.z;
    U[nc + 3] = st.rhoE;
    return U;
}

template<int nc>
Vec<nc> formFlux(const state<nc> &st, const vec &n, const gasinfo<nc> &gas) {
    Vec<nc> F;
    double rho = st.density();
    double vn = st.rhou.dot(n) / rho;
    double p = gas.pressure(st);
    F[0] = st.rhou.dot(n);
    for (int i = 1; i < nc; i++)
        F[i] = F[0] * st.rho[i] / rho;
    F[nc + 0] = st.rhou.x * vn + p * n.x;
    F[nc + 1] = st.rhou.y * vn + p * n.y;
    F[nc + 2] = st.rhou.z * vn + p * n.z;
    F[nc + 3] = (st.rhoE + p) * vn;
    return F;
}

template<int nc>
state<nc> formState(const Vec<nc> &U) {
    state<nc> st;
    st.rho[0] = U[0];
    for (int i = 1; i < nc; i++) {
        st.rho[i] = U[i];
        st.rho[0] -= U[i];
    }
    st.rhou = vec(U[nc], U[nc + 1], U[nc + 2]);
    st.rhoE = U[nc + 3];
    return st;
}

template<int nc>
Mat<nc> formOmega(const Vec<nc> &U, const vec &n, const gasinfo<nc> &gas) {
    state<nc> st = formState<nc>(U);
    double rho = st.density();
    vec v = st.rhou / rho;
    double vn = v.dot(n);
    double q2 = v.dot(v);
    double c = gas.sound_speed(st);
    double g = gas.gamma_factor(st);
    double b = 1 / (g - 1);
    double t[nc];
    st.fractions(t);

    Mat<nc> Om;
    Om.fill(0);
    for (int i = 1; i < nc; i++) {
        Om(i - 1, 0) = -U[i] / U[0];
        Om(i - 1, i) = 1;
    }

    Om(nc - 1, 0) = n.x > 0.5 ? vn * vn - q2/2 - c*c*b : -v.x;
    Om(nc    , 0) = n.y > 0.5 ? vn * vn - q2/2 - c*c*b : -v.y;
    Om(nc + 1, 0) = n.z > 0.5 ? vn * vn - q2/2 - c*c*b : -v.z;
    Om(nc - 1, nc    ) = n.x > 0.5 ? -vn : 1;
    Om(nc    , nc + 1) = n.y > 0.5 ? -vn : 1;
    Om(nc + 1, nc + 2) = n.z > 0.5 ? -vn : 1;
    Om(nc - 1, nc + 3) = n.x;
    Om(nc    , nc + 3) = n.y;
    Om(nc + 1, nc + 3) = n.z;

    Om(nc + 2, 0) = q2 / 2 - c * vn / g;
    Om(nc + 3, 0) = q2 / 2 + c * vn / g;

    double M = gas.molar_mass(st);
    for (int i = 1; i < nc; i++) {
        double dbdti = M * (
                (b - gas.beta[0]) / gas.molar[0] - (b - gas.beta[i]) / gas.molar[i]
            );

        Om(nc + 2, i) = -c * c / g * dbdti;
        Om(nc + 3, i) = -c * c / g * dbdti;
        Om(nc + 2, 0) += c * c / g * dbdti * U[i] / U[0];
        Om(nc + 3, 0) += c * c / g * dbdti * U[i] / U[0];
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

template<int nc>
Vec<nc> formLambda(const Vec<nc> &U, const vec &n, const gasinfo<nc> &gas) {
    state<nc> st = formState<nc>(U);
    double rho = st.density();
    vec v = st.rhou / rho;
    double vn = v.dot(n);
    double c = gas.sound_speed(st);

    Vec<nc> lam;
    for (int i = 0; i < nc + 2; i++)
        lam[i] = vn;

    lam[nc + 2] = vn + c;
    lam[nc + 3] = vn - c;

    return lam;
}

template<int nc>
void solver_flux<nc>::solve(const state<nc> &le, const state<nc> &ri, const vec &n, const gasinfo<nc> &gas) {
    Vec<nc> UL, UR, U, lam;
    Vec<nc> FL, FR, F;
    Mat<nc> Om, iOm;

    UL = formVector(le);
    UR = formVector(ri);

    FL = formFlux(le, n, gas);
    FR = formFlux(ri, n, gas);

    U = 0.5 * (UL + UR);

    Om = formOmega(U, n, gas);
    lam = formLambda(U, n, gas);
    iOm = Om.inverse();

    F = .5 * (FL + FR + iOm * lam.cwiseAbs().cwiseProduct(Om * (UL - UR)));

    fden[0] = F[0];
    for (int i = 1; i < nc; i++) {
        fden[i] = F[i];
        fden[0] -= F[i];
    }
    fmom.x = F[nc + 0];
    fmom.y = F[nc + 1];
    fmom.z = F[nc + 2];
    fener = F[nc + 3];
}

#include <Eigen/Core>

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

    Mat<nc> Om;
    Om.fill(0);
    for (int i = 1; i < nc; i++) {
        Om(i - 1, 0) = -U[i] / U[0];
        Om(i - 1, i) = 1;
    }

    Om(nc - 1, 0) = n.x > 0.5 ? -vn * vn + q2/2 + c*c*b : -v.x;
    Om(nc    , 0) = n.y > 0.5 ? -vn * vn + q2/2 + c*c*b : -v.y;
    Om(nc + 1, 0) = n.z > 0.5 ? -vn * vn + q2/2 + c*c*b : -v.z;
    Om(nc - 1, nc    ) = n.x > 0.5 ? vn : 1;
    Om(nc    , nc + 1) = n.y > 0.5 ? vn : 1;
    Om(nc + 1, nc + 2) = n.z > 0.5 ? vn : 1;
    Om(nc - 1, nc + 3) = -n.x;
    Om(nc    , nc + 3) = -n.y;
    Om(nc + 1, nc + 3) = -n.z;

    Om(nc + 2, 0) = q2 / 2 - c * vn * b;
    Om(nc + 3, 0) = q2 / 2 + c * vn * b;

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
Mat<nc> formiOmega(const Vec<nc> &U, const vec &n, const gasinfo<nc> &gas) {
    state<nc> st = formState<nc>(U);
    double rho = st.density();
    vec v = st.rhou / rho;
    double vn = v.dot(n);
    double q2 = v.dot(v);
    double c = gas.sound_speed(st);
    double g = gas.gamma_factor(st);
    double b = 1 / (g - 1);
    const double c2b = c*c*b;

    Mat<nc> iOm;

    Vec<nc> col, row;
    col[0] = 1;
    for (int i = 1; i < nc; i++)
        col[i] = U[i] / U[0];
    col[nc + 0] = v.x;
    col[nc + 1] = v.y;
    col[nc + 2] = v.z;
    col[nc + 3] = q2/2 + c2b;

    double M = gas.molar_mass(st);
    for (int i = 1; i < nc; i++) {
        double dbdti = M * (
                (b - gas.beta[0]) / gas.molar[0] - (b - gas.beta[i]) / gas.molar[i]
            );

        row[i - 1] = c * c / g * dbdti;
    }
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
    iOm = formiOmega(U, n, gas);
#ifndef NDEBUG
    const auto &Id = Om * iOm;
    const auto &Id2 = Mat<nc>::Identity();
    if (!Id.isApprox(Id2, 1e-6)) {
        state<nc> st = formState<nc>(U);
        std::cout << "beta = " << gas.beta_factor(st) << std::endl;
        std::cout << "dbdt1 = " << gas.dbetadtheta(st, 1) << std::endl;
        std::cout << "c = " << gas.sound_speed(st) << std::endl;
        std::cout << "t1 = " << U[1] / U[0] << std::endl;
        std::cout << "v = " << st.velocity() << std::endl;
        std::cout << "W = " << std::endl << Om << std::endl;
        std::cout << "iW = " << std::endl << iOm << std::endl;
        std::cout << "W . iW = " << std::endl << Id << std::endl;
        abort();
    }
#endif

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

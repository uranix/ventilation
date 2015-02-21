#ifndef __SLOPE_H__
#define __SLOPE_H__

#include <Eigen/Core>

#include "../include/state.h"

template<int nc> using Vec = Eigen::Matrix<double, nc + 4, 1>;
template<int nc> using Mat = Eigen::Matrix<double, nc + 4, nc + 4>;

template<int nc>
struct slope {
    Vec<nc> dW;

    double theta[nc];
    double dbdth[nc];
    double gamma;
    double c;
    vec v;
    dir::Direction dir;

    static Vec<nc> stateToVec(const state<nc> &st);
    static state<nc> vecToState(const Vec<nc> &U);
    slope();
    slope(const state<nc> &le, const state<nc> &ri, const dir::Direction dir, const gasinfo<nc> &gas);
    Mat<nc> Omega() const;
    Mat<nc> iOmega() const;
    Vec<nc> lambda() const;
};

#endif

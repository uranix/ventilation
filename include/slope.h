#ifndef __SLOPE_H__
#define __SLOPE_H__

#include <Eigen/Core>

#include "../include/state.h"

using Vec = Eigen::Matrix<double, nc + 4, 1>;
using Mat = Eigen::Matrix<double, nc + 4, nc + 4>;

struct slope {
    Vec dW;

    double theta[nc];
    double dbdth[nc];
    double gamma;
    double c;
    vec v;
    dir::Direction dir;

    static Vec stateToVec(const state &st);
    static Vec stateToFlux(const state &st, dir::Direction dir, const gasinfo &gas);
    static state vecToState(const Vec &U);
    slope();
    slope(const state &le, const state &ri, const dir::Direction dir, const gasinfo &gas);
    Mat Omega() const;
    Mat iOmega() const;
    Vec lambda() const;
    double lambda(size_t i) const;
};

#endif

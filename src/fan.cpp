#include "../include/fan.h"

namespace objects {

void fan::integrate_rhs(state &cell, const state &source, const double t, const double dt) {
    pipe::integrate_rhs(cell, source, t, dt);

    double s = this->surface * cell.velocity()(this->dir) / Qmax;
    double on = s > 1 ? exp(5 - 5*s) : 1;
    vec gradP = on * gradPmax;

    cell.rhou += dt * gradP;
    cell.e += dt * gradP.dot(cell.velocity());
}

}

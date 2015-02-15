#ifndef __PULSE_H__
#define __PULSE_H__

#include "object.h"

namespace objects {

template<int nc>
struct pulse : public object<nc> {
    const double drE, freq;
    const double pi;
    pulse(const vec &ll, const vec &ur, const std::string &id, const double drE, const double freq)
        : object<nc>(1, 1, 1, ll, ur, id), drE(drE), freq(freq), pi(4 * atan(1))
    { }

    virtual void compute_outer_flux(dir::Direction) override { }
    virtual void integrate_by(dir::Direction, const double, const double) override { }

    virtual double get_max_dt() const override { return .2 / freq; }

    virtual void integrate_rhs(state<nc> &cell, const state<nc> &source, const double t, const double dt) override {
        object<nc>::integrate_rhs(cell, source, t, dt);

        double A = sin(2 * pi * freq * (t + dt)) - sin(2 * pi * freq * t);
        this->ref(0, 0, 0).rhoE += drE * A;
    }
};

}

#endif

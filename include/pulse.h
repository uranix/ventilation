#ifndef __PULSE_H__
#define __PULSE_H__

#include "object.h"

namespace objects {

struct pulse : public object {
    const double de, freq;
    const double pi;
    pulse(const vec &ll, const vec &ur, const std::string &id, const double de, const double freq)
        : object(1, 1, 1, ll, ur, id), de(de), freq(freq), pi(4 * atan(1))
    { }

    virtual void compute_outer_flux(dir::Direction) override { }
    virtual void integrate_by(dir::Direction, const double, const double) override { }

    virtual double get_max_dt() const override { return .2 / freq; }

    virtual void integrate_rhs(state &cell, const state &source, const double t, const double dt) override {
        object::integrate_rhs(cell, source, t, dt);

        double A = sin(2 * pi * freq * (t + dt)) - sin(2 * pi * freq * t);
        this->ref(0, 0, 0).e += de * A;
    }
};

}

#endif

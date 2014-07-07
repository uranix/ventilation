#ifndef __SLOPED_STATE_H__
#define __SLOPED_STATE_H__

#include "state.h"

inline double touch(double value) {
    if (value < 0)
        return sin(value);
    else
        return cos(value);
}

template <int nc>
struct const_sloped_state {
    const state<nc> &avg;
    const state<nc> &sx;
    const state<nc> &sy;
    const state<nc> &sz;
    const_sloped_state(const state<nc> &avg, const state<nc> &sx, const state<nc> &sy, const state<nc> &sz)
        : avg(avg), sx(sx), sy(sy), sz(sz)
    {
    }
    const state<nc> &ce() const {
        return avg;
    }
    template<bool add>
    const state<nc> rec(const state<nc> &slope) {
        state<nc> ret = avg;
        if (add) {
            for (int i = 0; i < nc; i++)
                ret.rho[i] += slope.rho[i];
            ret.rhou += slope.rhou;
            ret.rhoE += slope.rhoE;
        } else {
            for (int i = 0; i < nc; i++)
                ret.rho[i] -= slope.rho[i];
            ret.rhou -= slope.rhou;
            ret.rhoE -= slope.rhoE;
        }
        return ret;
    }

    template<bool posdir>
    const state<nc> get(dir::Direction dir) {
        if (dir == dir::X)
            return rec<posdir>(sx);
        if (dir == dir::Y)
            return rec<posdir>(sy);
        return rec<posdir>(sz);
    }

    const state<nc> &slope(dir::Direction dir) const {
        if (dir == dir::X)
            return sx;
        if (dir == dir::Y)
            return sy;
        return sz;
    }
};

template <int nc>
struct sloped_state {
    state<nc> &avg;
    state<nc> &sx;
    state<nc> &sy;
    state<nc> &sz;
    sloped_state(state<nc> &avg, state<nc> &sx, state<nc> &sy, state<nc> &sz)
        : avg(avg), sx(sx), sy(sy), sz(sz)
    { }
    state<nc> &slope(dir::Direction dir) {
        if (dir == dir::X)
            return sx;
        if (dir == dir::Y)
            return sy;
        return sz;
    }
};

#endif

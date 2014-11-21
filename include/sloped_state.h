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
#if SECOND_ORDER
    const state<nc> &sx;
    const state<nc> &sy;
    const state<nc> &sz;
#endif
    const_sloped_state(const state<nc> &avg
#if SECOND_ORDER
            , const state<nc> &sx, const state<nc> &sy, const state<nc> &sz
#endif
        ) : avg(avg)
#if SECOND_ORDER
          , sx(sx), sy(sy), sz(sz)
#endif
    {
    }
    const state<nc> &ce() const {
        return avg;
    }

#if SECOND_ORDER
    template<bool posdir>
    const state<nc> rec(const state<nc> &slope) {
        state<nc> ret = avg;
        if (posdir) {
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
#endif

    template<bool posdir>
    const state<nc> get(dir::Direction dir) {
#if SECOND_ORDER
        return rec<posdir>(slope(dir));
#else
        (void)dir;
        return avg;
#endif
    }

#if SECOND_ORDER
    const state<nc> &slope(dir::Direction dir) const {
        if (dir == dir::X)
            return sx;
        if (dir == dir::Y)
            return sy;
        return sz;
    }
#endif
};

template <int nc>
struct sloped_state {
    state<nc> &avg;
#if SECOND_ORDER
    state<nc> &sx;
    state<nc> &sy;
    state<nc> &sz;
#endif
    sloped_state(state<nc> &avg
#if SECOND_ORDER
            , state<nc> &sx, state<nc> &sy, state<nc> &sz
#endif
            )
        : avg(avg)
#if SECOND_ORDER
          , sx(sx), sy(sy), sz(sz)
#endif
    { }
#if SECOND_ORDER
    state<nc> &slope(dir::Direction dir) {
        if (dir == dir::X)
            return sx;
        if (dir == dir::Y)
            return sy;
        return sz;
    }
#endif
};

#endif

#ifndef __ROOM_H__
#define __ROOM_H__

#include "object.h"

namespace objects {

template<int nc>
struct room : public object<nc> {
    room(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id)
        : object<nc>(nx, ny, nz, ll, ur, id)
    { }
    virtual ~room() { }
};

}

#endif

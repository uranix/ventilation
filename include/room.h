#ifndef __ROOM_H__
#define __ROOM_H__

#include "object.h"

namespace objects {

struct room : public object {
    room(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id)
        : object(nx, ny, nz, ll, ur, id)
    { }
    virtual ~room() { }
};

}

#endif

#ifndef __ROOM_H__
#define __ROOM_H__

#include "scene_object.h"

namespace objects {

template<int nc>
struct room : public scene_object<nc> {
	room(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id)
		: scene_object<nc>(nx, ny, nz, ll, ur, id)
	{
	}
	
	virtual ~room() {
	}
};

}

#endif

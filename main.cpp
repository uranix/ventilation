#include "room.h"
#include "pipe.h"
#include "solver.h"

#include <iostream>

#include <fenv.h>

int main() {

	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

	gasinfo<2> gas;

	gas.set(0, 29, 1.4); /* air */
	gas.set(1, 16, 1.33); /* methane */

	room<2> c(20, 20, 20, vec(0, 0, 0), vec(.3, .3, .3), "room", gas);
	pipe<2> c3(20, 0, vec(-1, .2, 0), vec(0, .26, .06), "pipe2", gas);

	double r1[2] = {1.20, 0};
	double r2[2] = {1.00, 0.03};

	for (int i = 0; i < c.nx; i++)
		for (int j = 0; j < c.ny; j++)
			for (int k = 0; k < c.nz; k++) {
				c(i, j, k).from_ruT(r2, vec(0, 0, 0), 300, gas);
			}
	
	for (int i = 0; i < c3.nx; i++)
		for (int j = 0; j < c3.ny; j++)
			for (int k = 0; k < c3.nz; k++) {
				c3(i, j, k).from_ruT(r2, vec(0, 0, 0), 300, gas);
			}

	connect(c, c3);

	double tmax = 0.1;
	double t = 0;
	
	int step = 0;

	std::vector<objects::scene_object<2> *> scene;
	scene.push_back(&c);
	scene.push_back(&c3);

	solver<2> solver(scene, 0.25);

	return 0;
}

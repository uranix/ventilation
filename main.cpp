#include <iostream>

#include "room.h"
#include "pipe.h"
#include "fan.h"
#include "atm.h"

#include "solver.h"

#include <fenv.h>

struct AtmValues : public functor<2> {
	const gasinfo<2> gas;
	AtmValues(const gasinfo<2> &gas) : gas(gas) { }
	void operator()(const vec &p, state<2> &st) const {
		std::vector<double> r = {1, 0};
		st.from_ruT(r, vec(0, 0, 0), 300, gas);
	}
};

struct RoomValues : public functor<2> {
	const gasinfo<2> gas;
	RoomValues(const gasinfo<2> &gas) : gas(gas) { }
	void operator()(const vec &p, state<2> &st) const {
		std::vector<double> r = {0.81875, 0};
		st.from_ruT(r, vec(0, 0, 0), 300, gas);
	}
};

int main() {
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

	gasinfo<2> gas;

	gas.set(0, 29, 1.4, 1.78e-5); /* air */
	gas.set(1, 16, 1.33, 1e-5); /* methane */

	objects::scene_object<2>::set_gas(gas);

	objects::room<2> room(20, 20, 20, vec(0, 0, 0), vec(1, 1, 1), "room");
	objects::fan<2> fan(10, 0, vec(-1, .45, .45), vec(0, .55, .55), "fan", 50000, 10);
	objects::pipe<2> pipe2(10, 1, vec(.9, 1, 0), vec(1, 1.5, .1), "pipe2");
	objects::atm<2> atm(vec(.8, 1.5, -.1), vec(1.1, 1.8, .2), "atm");
	objects::atm<2> atm2(vec(-1.2, .4, .4), vec(-1, .6, .6), "atm2");

	room.fill(RoomValues(gas));
	pipe2.fill(AtmValues(gas));
	atm.fill(AtmValues(gas));
	atm2.fill(AtmValues(gas));
	fan.fill(AtmValues(gas));

	std::vector<objects::scene_object<2> *> scene;
	scene.push_back(&fan);
	scene.push_back(&room);
	scene.push_back(&pipe2);
	scene.push_back(&atm);
	scene.push_back(&atm2);

	solver<2> solver(scene, 0.25);

	while (solver.time() < 10) {
		solver.compute_fluxes();
		double dt = solver.estimate_timestep();
		solver.integrate(dt);
		if ((solver.step() % 50) == 0) {
			std::cout << "t = " << solver.time() << ", dt = " << dt << ", step = " << solver.step() << std::endl;
			solver.save("vtks/");
		}
	}

	return 0;
}

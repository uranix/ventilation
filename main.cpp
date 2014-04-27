#include "room.h"

#include <iostream>

#include <fenv.h>

int main() {

	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

	gasinfo<2> gas;

	gas.set(0, 29, 1.4); /* air */
	gas.set(1, 16, 1.33); /* methane */

	room<2> c(20, 20, 20, vec(0, 0, 0), vec(.3, .3, .3), "room", gas);
	room<2> c2(20, 3, 3, vec(-1, 0, .2), vec(0, .06, .26), "pipe", gas);
//	room<2> c3(1, 20, 1, vec(.2, -1, .1), vec(.3, 0, .2), "pipe2", gas);
	room<2> c3(20, 3, 3, vec(-1, .2, 0), vec(0, .26, .06), "pipe2", gas);

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

	for (int i = 0; i < c2.nx; i++)
		for (int j = 0; j < c2.ny; j++)
			for (int k = 0; k < c2.nz; k++) {
				c2(i, j, k).from_ruT(r1, vec(10, 0, 0), 300, gas);
			}

	update_contacts(c, c2);
	update_contacts(c, c3);

	double tmax = 0.1;
	double t = 0;
	
	int step = 0;

	while (t < tmax) {
		c.compute_inner_fluxes();
		c.compute_outer_fluxes();
		
		c2.compute_inner_fluxes();
		c2.compute_outer_fluxes();

		c3.compute_inner_fluxes();
		c3.compute_outer_fluxes();

		double maxdt = c.get_min_h() / c.get_max_speed();
		double maxdt2 = c2.get_min_h() / c2.get_max_speed();
		double maxdt3 = c3.get_min_h() / c3.get_max_speed();
		
		double C = .25;
		double dt = C * std::min(std::min(maxdt, maxdt2), maxdt3);

		c.integrate(dt);
		c2.integrate(dt);
		c3.integrate(dt);

		if ((step % 100) == 0) {
			std::cout << "step = " << step << ", t = " << t << ", dt = " << dt << std::endl;
			c.save(step);
			c2.save(step);
			c3.save(step);
		}

		t += dt;
		step++;
	}

	return 0;
}

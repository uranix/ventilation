#ifndef __EULER_H__
#define __EULER_H__

#ifdef __CUDACC__
#define APICALL __device__ __host__
#else

#include <cmath>

#define APICALL
#endif

#ifndef __riemann_debug_stream

namespace riemann_solver {

enum Quiet { LOG };
template<typename T>
APICALL Quiet operator <<(Quiet o, T) { return o; }

}

#define __riemann_debug_stream LOG

#endif

namespace riemann_solver {

/** Riemann solver for 1D Euler equations. Adiabatic index on the left and on the right can be different */
template<typename real>
class euler_1d {
	real density[2];
	real velocity[2];
	real pressure[2];
	real gamma[2];

	real chi[2]; /* (g - 1) / 2 / g */
	real sound[2];

	real U, P;
	real R1, E1, R2, E2;

	real W1, W1s, W2s, W2;

	const real eps;
	const int max_iters;

	struct multivalue {
		real val, der;
		APICALL multivalue(real v, real d) : val(v), der(d) { }
	};
	
	template<int n>
	APICALL multivalue f(real P) const {
		real pi = P / pressure[n];

		if (pi == 0) {
			return multivalue(2 * sound[n] / (1 - gamma[n]), 0.);
		}

		real beta = sqrt((gamma[n] + 1) / (2 * gamma[n]) * pi + chi[n]);
		if (pi >= 1) {
			real b2 = beta * beta;
			real b3 = beta * b2;
			return multivalue(
				(P - pressure[n]) / density[n] / sound[n] / beta,
				((gamma[n] + 1) * pi + 3 * gamma[n] - 1) / 
					(4 * density[n] * sound[n] * gamma[n] * b3)
			);
		} else {
			real pix = pow(pi, chi[n]);
			return multivalue(
				2 * sound[n] / (gamma[n] - 1) * (pix - 1),
				sound[n] / gamma[n] / P * pix
			);
		}
	}

	APICALL multivalue F(real P) const {
		multivalue v1(f<0>(P));
		multivalue v2(f<1>(P));

		return multivalue(
			v1.val + v2.val - velocity[0] + velocity[1], 
			v1.der + v2.der);
	}

/* {{{ Initial condition approximation */
	APICALL real acoustic_approximation() const {
		real a1 = density[0] * sound[0];
		real a2 = density[1] * sound[1];

		real p1 = pressure[0];
		real p2 = pressure[1];

		real v = velocity[0] - velocity[1];

		return (p1 * a2 + p2 * a1 + v * a1 * a2) / (a1 + a2);
	}

	APICALL real rare_gas_approximation() const {
		real a1 = 2 * sound[0] / (gamma[0] - 1);
		real a2 = 2 * sound[1] / (gamma[1] - 1);

		real u1 = velocity[0];
		real u2 = velocity[1];

		real A = a1 / pow(pressure[0], chi[0]);
		real B = a2 / pow(pressure[1], chi[1]);

		real ichiavg = 1 / sqrt(chi[0] * chi[1]);

		return pow((u1 - u2 + a1 + a2) / (A + B), ichiavg);
	}
	
	APICALL real approximate_pressure() const {
		real P0 = acoustic_approximation();

		if (P0 <= 0)
			P0 = rare_gas_approximation();
		
		while (F(P0).val > 0)
			P0 /= 2;

		return P0;
	} /* }}} */

	APICALL bool solve_vacuum_case() {
		if (F(0).val < 0) 
			return false;

		R1 = R2 = E1 = E2 = P = 0;

		W1 = velocity[0] - sound[0];
		W1s = velocity[0] + 2 / (gamma[0] - 1) * sound[0];

		W2 = velocity[1] + sound[1];
		W2s = velocity[1] - 2 / (gamma[1] - 1) * sound[1];

		U = (W1s + W2s) / 2;

		__riemann_debug_stream << "Vacuum case, no iterations needed\n";

		return true;
	}

	APICALL static real max(real a, real b) {
		if (a < b)
			return a;
		return b;
	}

	APICALL bool newton_refine(real &P) {
		real pmax = max(pressure[0], pressure[1]);
		real dp;
		real safe_factor = static_cast<real>(0.99);
		int iters = 0;
		bool converged = false;
		do {
			multivalue v(F(P));
			dp = v.val / v.der;
			if (P - dp <= 0)
				dp = safe_factor * P;
			P -= dp;
			iters++;
			converged = fabs(dp) < eps * max(P, pmax);
		} while ((!converged) && (iters < max_iters));

		if (converged)
			__riemann_debug_stream << "Converged in " << iters << " Newton iterations\n";
		else
			__riemann_debug_stream << "Did not converge, done " << iters 
				<< " Newton iterations, but got only " << (dp / max(P, pmax)) 
				<< " relative error when " << eps << " requested\n";

		return converged;
	}

	APICALL void compute_waves() {
		if (P > pressure[0]) {
			lwave = SHOCK;

			real m1 = sqrt(density[0] / 2 * (P * (gamma[0] + 1) + pressure[0] * (gamma[0] - 1)));
			W1 = W1s = velocity[0] - m1 / density[0];
			R1 = density[0] * m1 / (m1 - (velocity[0] - U) * density[0]);
		} else {
			lwave = RAREFACTION;

			W1 = velocity[0] - sound[0];
			W1s = U - sound[0] + (gamma[0] - 1) / 2 * (U - velocity[0]);
			R1 = gamma[0] * P / (U - W1s) / (U - W1s);
		}
		E1 = P / R1 / (gamma[0] - 1);

		if (P > pressure[1]) {
			rwave = SHOCK;

			real m2 = sqrt(density[1] / 2 * (P * (gamma[1] + 1) + pressure[1] * (gamma[1] - 1)));

			W2 = W2s = velocity[1] + m2 / density[1];
			R2 = density[1] * m2 / (m2 + (velocity[1] - U) * density[1]);
		} else {
			rwave = RAREFACTION;

			W2 = velocity[1] + sound[1];
			W2s = U + sound[1] + (gamma[1] - 1) / 2 * (U - velocity[1]);
			R2 = gamma[1] * P / (W2s - U) / (W2s - U);
		}
		E2 = P / R2 / (gamma[1] - 1);
	}

	APICALL bool solve() {
		if (solve_vacuum_case()) {
			lwave = rwave = VACUUM;
			return true;
		}

		P = approximate_pressure();
		bool converged = newton_refine(P);

		U = (velocity[0] + velocity[1] - f<0>(P).val + f<1>(P).val) / 2;

		compute_waves();

		return converged;
	}
	
	APICALL static inline void safe_put(real A[], int i, real v) {
		if (A)
			A[i] = v;
	}

public:

	enum WaveType {
		VACUUM = 0,
		RAREFACTION = 1,
		SHOCK = 2
	};

	struct WaveFan {
		real W1, W1s, U, W2, W2s;
		APICALL WaveFan(real W1, real W1s, real U, real W2, real W2s) : W1(W1), W1s(W1s), U(U), W2(W2), W2s(W2s) { }
	};

private:
	WaveType lwave, rwave;
	bool solved_ok;
public:

	/** Solve problem with for specific values. 
	 * \param eps Relative error tolerance 
	 * \param max_iters Maximum number of Newton iteration to be performed. 
	 * If this number is exceeded iterations are stopped and solution is considered
	 * to be not good(), i.e. have error > eps.
	 */
	APICALL euler_1d(
		real density1, real velocity1, real pressure1, real gamma1,
		real density2, real velocity2, real pressure2, real gamma2,
		real eps = static_cast<real>(1e-12),
		int max_iters = 10
		)
	:
		eps(eps), max_iters(max_iters)
	{
		density[0] = density1;
		density[1] = density2;

		velocity[0] = velocity1;
		velocity[1] = velocity2;

		pressure[0] = pressure1;
		pressure[1] = pressure2;

		gamma[0] = gamma1;
		gamma[1] = gamma2;

		for (int n = 0; n < 2; n++) {
			chi[n] = (gamma[n] - 1) / (2 * gamma[n]);
			sound[n] = sqrt(gamma[n] * pressure[n] / density[n]);
		}

		solved_ok = solve();
	}

	/** Returns whether Newton iterations have converged to given precision or not */
	APICALL bool good() const {
		return solved_ok;
	}

	/** Returns left wave type - shock or rarification wave. 
	 *
	 * Special value ::VACUUM is used for case with vacuum between waves */
	APICALL WaveType left_wave_type() const {
		return lwave;
	}

	/** Returns right wave type - shock or rarification wave. 
	 *
	 * Special value ::VACUUM is used for case with vacuum between waves */
	APICALL WaveType right_wave_type() const {
		return rwave;
	}
	
	/** Return waves speed 
	 *
	 * Each rarification wave is represented as pair (Wn, Wns), 
	 * each shock wave has Wn = Wns. Contact discontinuity has velocity U */
	APICALL const WaveFan get_wave_fan() const {
		return WaveFan(W1, W1s, U, W2s, W2);
	}

	/** Return maximum absolute speed
	 *
	 * Can be used to estimate Courant's number */
	APICALL const real max_speed() const {
		return max(fabs(W1), fabs(W2));
	}

	/** Get solution as arrays of density, velocity, specific energy and pressure
	* 
	* \param N Number of elements in each array
	* \param xi Array with x/t values (x = 0 is the initial discontinuity).
	* \param dens Optional output array used for density values
	* \param vel Optional output array used for velocity values
	* \param specen Optional output array used for specific energy values
	* \param press Optional output array used for pressure values
	* \param incvel Velocity array storage increment. Intended for use in multidimensional solvers
	**/
	APICALL void get(const int N, const real xi[], real dens[], real vel[], real specen[], real press[], 
		const int incvel = 1) const 
	{
		for (int i = 0, iv = 0; i < N; i++, iv += incvel) {
			bool leftmost = xi[i] <= W1;
			bool rightmost = xi[i] >= W2;
			bool lefthalf = xi[i] < U;
			bool mid = xi[i] >= W1s && xi[i] <= W2s;
			bool leftmid = mid && lefthalf;
			bool rightmid = mid && !lefthalf;

			if (leftmost) {
				const int j = 0;
				safe_put(dens, i, density[j]);
				safe_put(vel, iv, velocity[j]);
				safe_put(specen, i, pressure[j] / density[j] / (gamma[j] - 1));
				safe_put(press, i, pressure[j]);
				continue;
			}
			if (rightmost) {
				const int j = 1;
				safe_put(dens, i, density[j]);
				safe_put(vel, iv, velocity[j]);
				safe_put(specen, i, pressure[j] / density[j] / (gamma[j] - 1));
				safe_put(press, i, pressure[j]);
				continue;
			}
			if (leftmid) {
				safe_put(dens, i, R1);
				safe_put(vel, iv, U);
				safe_put(specen, i, E1);
				safe_put(press, i, P);
				continue;
			}
			if (rightmid) {
				safe_put(dens, i, R2);
				safe_put(vel, iv, U);
				safe_put(specen, i, E2);
				safe_put(press, i, P);
				continue;
			}
			if (lefthalf) {
				/* W1 < xi < W1s */
				const int j = 0;
				real cs = 2 / (gamma[j] + 1) * sound[j] + (gamma[j] - 1) / (gamma[j] + 1) * (velocity[j] - xi[i]);
				real rs = density[j] * pow(cs / sound[j], 2 / (gamma[j] - 1));
				real ps = pressure[j] * pow(cs / sound[j], 2 * gamma[j] / (gamma[j] - 1));
				real es = ps / rs / (gamma[j] - 1);
				safe_put(dens, i, rs);
				safe_put(vel, iv, xi[i] + cs);
				safe_put(specen, i, es);
				safe_put(press, i, ps);
			} else {
				/* W2s < xi < W2 */
				const int j = 1;
				real cs = 2 / (gamma[j] + 1) * sound[j] - (gamma[j] - 1) / (gamma[j] + 1) * (velocity[j] - xi[i]);
				real rs = density[j] * pow(cs / sound[j], 2 / (gamma[j] - 1));
				real ps = pressure[j] * pow(cs / sound[j], 2 * gamma[j] / (gamma[j] - 1));
				real es = ps / rs / (gamma[j] - 1);
				safe_put(dens, i, rs);
				safe_put(vel, iv, xi[i] - cs);
				safe_put(specen, i, es);
				safe_put(press, i, ps);
			}
		}
	}
};


template<typename real, int dims>
class euler {
public:

	struct vec {
		real v[dims];
		APICALL real &operator[](int i) { return v[i]; }
		APICALL const real &operator[](int i) const { return v[i]; }
		APICALL real *data() { return &v[0]; }
		APICALL const real *data() const { return &v[0]; }
	};

private:
	vec velocity1;
	vec velocity2;
	vec dir;

	euler_1d<real> onedim_solver;

	APICALL static real dot(vec x, vec y) {
		real ret = 0;
		for (int i = 0; i < dims; i++)
			ret += x[i] * y[i];
		return ret;
	}

	APICALL static vec normalized(vec x) {
		real norm = sqrt(dot(x, x));
		vec ret;
		for (int i = 0; i < dims; i++)
			ret[i] = x[i] / norm;
		return ret;
	}

public:

	APICALL euler(
		real density1, vec velocity1, real pressure1, real gamma1,
		real density2, vec velocity2, real pressure2, real gamma2,
		vec direction,
		real eps = static_cast<real>(1e-12), int max_iters = 10)
	:
		velocity1(velocity1), velocity2(velocity2), dir(normalized(direction)),
		onedim_solver(
			density1, dot(velocity1, dir), pressure1, gamma1,
			density2, dot(velocity2, dir), pressure2, gamma2,
			eps, max_iters)
	{
	}

	APICALL const euler_1d<real> &get_1d_solver() const {
		return onedim_solver;
	}
	
	APICALL void get(const int N, const real xi[], real dens[], vec vel[], 
		real specen[], real press[]) const 
	{
		real U = onedim_solver.get_wave_fan().U;
		onedim_solver.get(N, xi, dens, reinterpret_cast<real *>(vel), specen, press, sizeof(vec) / sizeof(real));

		real v1n = dot(velocity1, dir);
		real v2n = dot(velocity2, dir);

		for (int i = 0; i < N; i++) {
			real vn = vel[i][0];
			if (xi[i] < U)
				for (int k = 0; k < dims; k++)
					vel[i][k] = (vn - v1n) * dir[k] + velocity1[k];
			else
				for (int k = 0; k < dims; k++)
					vel[i][k] = (vn - v2n) * dir[k] + velocity2[k];
		}
	}
};

}

#endif

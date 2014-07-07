#include "box.h"

#include <cmath>
#include <iostream>

double intersect_rect(dir::Direction d, const vec &all, const vec &aur,
    const vec &bll, const vec &bur, double tol)
{
    vec ins;

    ins.x = std::min(aur.x, bur.x) - std::max(all.x, bll.x);
    ins.y = std::min(aur.y, bur.y) - std::max(all.y, bll.y);
    ins.z = std::min(aur.z, bur.z) - std::max(all.z, bll.z);

    double S = 1;
	for (auto dd : dir::DIRECTIONS) {
		if (dd == d)
			continue;
        if (ins(dd) < tol)
            return -1;
        else
            S *= ins(dd);
	}
    return S;
}

double try_side(dir::Direction d, int s, const box &a, const box &b, double tol) {
    if (a.closed[d][s] || b.closed[d][1-s])
        return -1;

    double ax = (s == 0) ? a.ll(d) : a.ur(d);
    double bx = (s == 1) ? b.ll(d) : b.ur(d);

    if (fabs(ax - bx) > tol)
        return -1;

    return intersect_rect(d, a.ll, a.ur, b.ll, b.ur, tol);
}

void connect(box &a, box &b) {
    double unit_length = std::min((a.ur - a.ll).norm(), (b.ur - b.ll).norm());
    double tol = 1e-6;

    tol *= unit_length;

    bool found = false;

	for (auto d : dir::DIRECTIONS)
        for (int s = 0; s < 2; s++) {
            bool check = try_side(d, s, a, b, tol) > 0;
            if (check && found) {
                std::cerr << "Multiple contacts! This should not have happened." << std::endl;
                return;
            }
            found |= check;

            if (!check)
                continue;

            std::cout << "Found contact in " << dir::to_char(d)
				<< " direction between " << (s ? a.id : b.id) << " -> " << (s ? b.id : a.id) << std::endl;

            double Sa = a.h.x * a.h.y * a.h.z / a.h(d);
            double Sb = b.h.x * b.h.y * b.h.z / b.h(d);

			int alo[3];
			int ahi[3];
			int blo[3];
			int bhi[3];
			for (auto dd : dir::DIRECTIONS) {
				alo[dd] = (d != dd) ? 0 : (s * (a.n(dd) - 1));
				ahi[dd] = (d != dd) ? a.n(dd) : (s * (a.n(dd) - 1) + 1);
				blo[dd] = (d != dd) ? 0 : ((1 - s) * (b.n(dd) - 1));
				bhi[dd] = (d != dd) ? b.n(dd) : ((1 - s) * (b.n(dd) - 1) + 1);
			}

#define LOOP(x, dir) for (int x ## dir = x ## lo[dir]; x ## dir < x ## hi[dir]; x ## dir ++)

            LOOP(a, 0)
            LOOP(a, 1)
            LOOP(a, 2)
            LOOP(b, 0)
            LOOP(b, 1)
            LOOP(b, 2)
            {
                vec all(a.ll);
                vec bll(b.ll);
                all += a.h * vec(a0, a1, a2);
                bll += b.h * vec(b0, b1, b2);
                vec aur = all + a.h;
                vec bur = bll + b.h;

                double S = intersect_rect(d, all, aur, bll, bur, tol);

                if (S > 0) {
                    a.side(d, s, a0, a1, a2).push_back(
                        box::Contact(&b, S / Sa, b0, b1, b2));
                    b.side(d, 1 - s, b0, b1, b2).push_back(
                        box::Contact(&a, S / Sb, a0, a1, a2));
                }
            }
#undef LOOP
        }
}

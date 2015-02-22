#include "box.h"

#include <cmath>
#include <iostream>
#include <cstdlib>

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

double try_side(dir::Direction d, dir::Side s, const box &a, const box &b, double tol) {
    if (a.closed[d][s] || b.closed[d][flip(s)])
        return -1;

    double ax = (s == dir::BEG) ? a.ll(d) : a.ur(d);
    double bx = (s == dir::END) ? b.ll(d) : b.ur(d);

    if (fabs(ax - bx) > tol)
        return -1;

    return intersect_rect(d, a.ll, a.ur, b.ll, b.ur, tol);
}

void connect(box &a, box &b) {
    double unit_length = std::min((a.ur - a.ll).norm(), (b.ur - b.ll).norm());
    double tol = 1e-6;

    tol *= unit_length;

    bool found = false;
    bool ok = false;

    for (auto d : dir::DIRECTIONS)
        for (auto s : dir::SIDES) {
            bool check = try_side(d, s, a, b, tol) > 0;
            if (check && found) {
                std::cerr << "Multiple contacts! This should not have happened." << std::endl;
                abort();
            }
            found |= check;

            if (!check)
                continue;

            std::cout << "Found contact in "
                << dir::to_char(d)
                << " direction between "
                << (s == dir::BEG ? a.id : b.id)
                << " -> "
                << (s == dir::END ? a.id : b.id)
                << std::endl;

            double Sa = a.h.x * a.h.y * a.h.z / a.h(d);
            double Sb = b.h.x * b.h.y * b.h.z / b.h(d);

            int alo[dir::DIR_END];
            int ahi[dir::DIR_END];
            int blo[dir::DIR_END];
            int bhi[dir::DIR_END];
            for (auto dd : dir::DIRECTIONS) {
                if (d == dd) {
                    alo[dd] = (s == dir::BEG) ? 0 : (a.n(dd) - 1);
                    ahi[dd] = (s == dir::BEG) ? 1 : a.n(dd);
                    blo[dd] = (s == dir::END) ? 0 : (b.n(dd) - 1);
                    bhi[dd] = (s == dir::END) ? 1 : b.n(dd);
                } else {
                    alo[dd] = 0;
                    ahi[dd] = a.n(dd);
                    blo[dd] = 0;
                    bhi[dd] = b.n(dd);
                }
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
                    b.side(d, flip(s), b0, b1, b2).push_back(
                        box::Contact(&a, S / Sb, a0, a1, a2));
                }
            }
#undef LOOP
            ok = true;
        }
    if (!ok) {
        std::cerr << "No contact between " << a.id << " and " << b.id << " found!" << std::endl;
        abort();
    }
}

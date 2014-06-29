#include "box.h"

#include <cmath>
#include <iostream>

double vc(const vec &v, int i) {
    if (i == 0)
        return v.x;
    if (i == 1)
        return v.y;
    return v.z;
}

double intersect_rect(int d, const vec &all, const vec &aur, 
    const vec &bll, const vec &bur, double tol) 
{
    vec ins;

    ins.x = std::min(aur.x, bur.x) - std::max(all.x, bll.x);
    ins.y = std::min(aur.y, bur.y) - std::max(all.y, bll.y);
    ins.z = std::min(aur.z, bur.z) - std::max(all.z, bll.z);

    double S = 1;
    if (d != 0) {
        if (ins.x < tol)
            return -1;
        else
            S *= ins.x;
    }
    if (d != 1) {
        if (ins.y < tol)
            return -1;
        else
            S *= ins.y;
    }
    if (d != 2) {
        if (ins.z < tol)
            return -1;
        else
            S *= ins.z;
    }

    return S;
}

double try_side(int d, int s, const box &a, const box &b, double tol) {
    if (a.closed[d][s] || b.closed[d][1-s])
        return -1;
    
    double ax = (s == 0) ? vc(a.ll, d) : vc(a.ur, d);
    double bx = (s == 1) ? vc(b.ll, d) : vc(b.ur, d);

    if (fabs(ax - bx) > tol)
        return -1;
    
    return intersect_rect(d, a.ll, a.ur, b.ll, b.ur, tol);
}

void connect(box &a, box &b) {
    double unit_length = std::min((a.ur - a.ll).norm(), (b.ur - b.ll).norm());
    double tol = 1e-6;

    tol *= unit_length;

    bool found = false;

    for (int d = 0; d < 3; d++)
        for (int s = 0; s < 2; s++) {
            bool check = try_side(d, s, a, b, tol) > 0;
            if (check && found) {
                std::cerr << "Multiple contacts! This should not have happened." << std::endl;
                return;
            }
            found |= check;

            if (!check)
                continue;

            char cdir = d == 0 ? 'X' : (d == 1 ? 'Y' : 'Z');
            std::cout << "Found contact in " << cdir << " direction between " << (s ? a.id : b.id) <<
                " -> " << (s ? b.id : a.id) << std::endl;

            double Sa = a.h.x * a.h.y * a.h.z / vc(a.h, d);
            double Sb = b.h.x * b.h.y * b.h.z / vc(b.h, d);

            int i1lo = (d != 0) ? 0 : (s * (a.nx - 1));
            int j1lo = (d != 1) ? 0 : (s * (a.ny - 1));
            int k1lo = (d != 2) ? 0 : (s * (a.nz - 1));
            int i1hi = (d != 0) ? a.nx : (s * (a.nx - 1) + 1);
            int j1hi = (d != 1) ? a.ny : (s * (a.ny - 1) + 1);
            int k1hi = (d != 2) ? a.nz : (s * (a.nz - 1) + 1);
            
            int i2lo = (d != 0) ? 0 : ((1 - s) * (b.nx - 1));
            int j2lo = (d != 1) ? 0 : ((1 - s) * (b.ny - 1));
            int k2lo = (d != 2) ? 0 : ((1 - s) * (b.nz - 1));
            int i2hi = (d != 0) ? b.nx : ((1 - s) * (b.nx - 1) + 1);
            int j2hi = (d != 1) ? b.ny : ((1 - s) * (b.ny - 1) + 1);
            int k2hi = (d != 2) ? b.nz : ((1 - s) * (b.nz - 1) + 1);

#define LOOP(x) for (int x = x ## lo; x < x ## hi; x++)

            LOOP(i1)
            LOOP(j1)
            LOOP(k1)
            LOOP(i2)
            LOOP(j2)
            LOOP(k2)
            {
                vec all(a.ll);
                vec bll(b.ll);
                all += a.h * vec(i1, j1, k1);
                bll += b.h * vec(i2, j2, k2);
                vec aur = all + a.h;
                vec bur = bll + b.h;

                double S = intersect_rect(d, all, aur, bll, bur, tol);

                if (S > 0) {
                    a.side(d, s, i1, j1, k1).push_back(
                        box::Contact(&b, S / Sa, i2, j2, k2));
                    b.side(d, 1 - s, i2, j2, k2).push_back(
                        box::Contact(&a, S / Sb, i1, j1, k1));
                }
            }
#undef LOOP
        }
}

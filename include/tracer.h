#ifndef __TRACER_H__
#define __TRACER_H__

#include "vec.h"

#include <vector>
#include <string>
#include <fstream>
#include <functional>

typedef std::pair<double, vec> PathSpec;

class tracer {
    std::vector<PathSpec> path;
    std::string name;
public:
    tracer(const std::vector<vec> &chain, const std::string &name, double step) : name(name) {
        path.push_back(PathSpec(0, chain.front()));
        for (size_t i = 1; i < chain.size(); i++) {
            int c = ceil((chain[i] - chain[i-1]).norm() / step);
            vec h = (chain[i] - chain[i-1]) / c;
            double step = h.norm();
            double start = path.back().first;
            for (int j = 1; j <= c; j++)
                path.push_back(PathSpec(start + j * step, chain[i-1] + j * h));
        }
    }
    template<int nc>
    void walk(const std::string &prefix, const int step, const gasinfo<nc> gas, std::function<state<nc>(vec)> locator) {
        std::string fn(prefix + name);
        fn += ".";
        fn += std::to_string(step);
        fn += ".csv";

        std::fstream f(fn.c_str(), std::ios::out);

        if (!f) {
            std::cerr << "Could not save trace " << fn << std::endl;
            return;
        }

        f << "len,";
        for (int i = 0; i < nc; i++)
            f << "rho" << i << ",";
        f << "rho,gamma,ux,uy,uz,v,p,eps,T" << std::endl;
        f.precision(15);

        size_t ps = path.size();
        for (size_t i = 0; i < ps; i++) {
            vec tang;
            if (i == 0)
                tang = path[1].second - path[0].second;
            else
                tang = path[i].second - path[i-1].second;

            tang *= 1. / tang.norm();

            const state<nc> &st = locator(path[i].second);

            f << path[i].first << ",";
            for (int j = 0; j < nc; j++)
                f << st.rho[j] << ",";
            f << st.density() << "," << gas.gamma_factor(st) << ",";
            vec v = st.velocity();
            f << v.x << "," << v.y << "," << v.z << "," << v.dot(tang) << ",";
            f << gas.pressure(st) << "," << st.specific_energy() << ","
                << gas.temperature(st) << std::endl;
        }
    }
};

#endif

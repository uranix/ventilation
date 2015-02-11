#include "../include/solver.h"

template class solver<NC>;

template<int nc>
solver<nc>::solver(
        const std::vector<objects::object<nc> *> &scene,
        const std::vector<tracer *> &tracers, const double cou)
: scene(scene), tracers(tracers), cou(cou)
{
#if FP_TRAP
    feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif
    for (auto p : scene)
        p->set_solver(this);
    t = 0;
    _step = 0;
}

template<int nc>
void solver<nc>::save(const std::string &prefix) {
    for (auto p : scene)
        p->save(prefix, _step);
    for (auto t : tracers) {
        t->template walk<nc>(prefix, _step, gas(), [this] (vec p) -> const state<nc> {
            for (auto o : scene) {
                if (o->has_point(p))
                    return o->state_at(p);
            }
            state<nc> wrong;
            wrong.rho[0] = 1;
            return wrong;
        });
    }
}

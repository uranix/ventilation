
    MAYBECONST state<nc> &operator()(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);
        return _states[i + (j + k * ny) * nx];
    }

    MAYBECONST state<nc> &operator()(dir::Direction dir, int i) MAYBECONST {
        if (dir == dir::X)
            return (*this)(i, 0, 0);
        if (dir == dir::Y)
            return (*this)(0, i, 0);
        return (*this)(0, 0, i);
    }

    MAYBECONST state<nc> &source(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);
        return _sources[i + (j + k * ny) * nx];
    }

    MAYBECONST flux<nc> &x_flux(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i <= nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);
        return _fluxes[0][i + (j + k * ny) * (nx + 1)];
    }

    MAYBECONST flux<nc> &y_flux(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j <= ny);
        assert(k >= 0 && k < nz);
        return _fluxes[1][i + (j + k * (ny + 1)) * nx];
    }

    MAYBECONST flux<nc> &z_flux(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k <= nz);
        return _fluxes[2][i + (j + k * ny) * nx];
    }

    MAYBECONST flux<nc> &flux_dir(const dir::Direction dir, int i) MAYBECONST {
        if (dir == dir::X)
            return x_flux(i, 0, 0);
        if (dir == dir::Y)
            return y_flux(0, i, 0);
        return z_flux(0, 0, i);
    }

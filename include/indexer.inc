// vim:ft=cpp et sts=4 ts=4 sw=4 si

    MAYBECONST state &source(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);
        return _sources[i + (j + k * ny) * nx];
    }

    MAYBECONST flux &x_flux(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i <= nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);
        return _fluxes[0][i + (j + k * ny) * (nx + 1)];
    }

    MAYBECONST flux &y_flux(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j <= ny);
        assert(k >= 0 && k < nz);
        return _fluxes[1][i + (j + k * (ny + 1)) * nx];
    }

    MAYBECONST flux &z_flux(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k <= nz);
        return _fluxes[2][i + (j + k * ny) * nx];
    }

    MAYBECONST slope &x_slope(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i <= nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);
        return _slopes[0][i + (j + k * ny) * (nx + 1)];
    }

    MAYBECONST slope &y_slope(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j <= ny);
        assert(k >= 0 && k < nz);
        return _slopes[1][i + (j + k * (ny + 1)) * nx];
    }

    MAYBECONST slope &z_slope(int i, int j, int k) MAYBECONST {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k <= nz);
        return _slopes[2][i + (j + k * ny) * nx];
    }

    MAYBECONST flux &flux_by(const dir::Direction dir, int i, int j, int k) MAYBECONST {
        if (dir == dir::X)
            return x_flux(i, j, k);
        if (dir == dir::Y)
            return y_flux(i, j, k);
        return z_flux(i, j, k);
    }

    MAYBECONST slope &slope_by(const dir::Direction dir, int i, int j, int k) MAYBECONST {
        if (dir == dir::X)
            return x_slope(i, j, k);
        if (dir == dir::Y)
            return y_slope(i, j, k);
        return z_slope(i, j, k);
    }

    MAYBECONST flux &flux_by(const dir::Direction dir, int i) MAYBECONST {
        if (dir == dir::X)
            return x_flux(i, 0, 0);
        if (dir == dir::Y)
            return y_flux(0, i, 0);
        return z_flux(0, 0, i);
    }

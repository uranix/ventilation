template<int nc>
template<typename T>
void scene_object<nc>::put(std::fstream &f, T value) const {
    char buf[sizeof(T)];
    *reinterpret_cast<T *>(&buf[0]) = value;
    for (int s = 0; s < sizeof(T) / 2; s++) {
        std::swap(buf[s], buf[sizeof(T) - 1 - s]);
    }
    f.write(buf, sizeof(T));
}

template<int nc>
void scene_object<nc>::save(const std::string &prefix, const int step) const {
    std::string fn(prefix + id);
    fn += ".";
    fn += std::to_string(step);
    fn += ".vtk";

    std::fstream f(fn.c_str(), std::ios::out | std::ios::binary);

    if (!f) {
        std::cerr << "Could not save file " << fn << std::endl;
        return;
    }

    f << "# vtk DataFile Version 3.0\n";
    f << "Block dump\n";
    f << "BINARY\n";
    f << "DATASET RECTILINEAR_GRID\n";
    f << "DIMENSIONS " << nx + 1 << " " << ny + 1 << " " << nz + 1;
    f << "\nX_COORDINATES " << nx + 1 << " float\n";
    for (int i = 0; i <= nx; i++)
        put<float>(f, ll.x + i * h.x);
    f << "\nY_COORDINATES " << ny + 1 << " float\n";
    for (int j = 0; j <= ny; j++)
        put<float>(f, ll.y + j * h.y);
    f << "\nZ_COORDINATES " << nz + 1 << " float\n";
    for (int k = 0; k <= nz; k++)
        put<float>(f, ll.z + k * h.z);
    f << "\nCELL_DATA " << nx*ny*nz;
    for (int ic = 0; ic < nc; ic++) {
        f << "\nSCALARS rho" << ic << " float\nLOOKUP_TABLE default\n";
        for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                    put<float>(f, (*this)(i, j, k).ce().rho[ic]);
    }
    f << "\nVECTORS v float\n";
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                const vec &v = (*this)(i, j, k).ce().velocity();
                put<float>(f, v.x);
                put<float>(f, v.y);
                put<float>(f, v.z);
            }
    f << "\nSCALARS p float\nLOOKUP_TABLE default\n";
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                put<float>(f, gas().pressure((*this)(i, j, k).ce()));
    f << "\nSCALARS T float\nLOOKUP_TABLE default\n";
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                put<float>(f, gas().temperature((*this)(i, j, k).ce()));
    f << "\nSCALARS eps float\nLOOKUP_TABLE default\n";
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                put<float>(f, gas().specific_energy((*this)(i, j, k).ce()));
    f << std::endl;
    f.close();
}

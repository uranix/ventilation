#include "../include/object.h"
#include <fstream>

namespace objects {

template struct object<NC>;

template<typename T>
void put(std::fstream &f, T value) {
    union {
        char buf[sizeof(T)];
        T val;
    } helper;
    helper.val = value;
    std::reverse(helper.buf, helper.buf + sizeof(T));
    f.write(helper.buf, sizeof(T));
}

template<int nc>
void object<nc>::save(const std::string &prefix, const int step) const {
    #include "GitVersion.h"
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
    f << "Block dump. " << VERSION << "\n";
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
                    put<float>(f, val(i, j, k).rho[ic]);
    }
    f << "\nSCALARS rho float\nLOOKUP_TABLE default\n";
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                put<float>(f, val(i, j, k).density());
    f << "\nSCALARS gamma float\nLOOKUP_TABLE default\n";
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                put<float>(f, gas().gamma_factor(val(i, j, k)));
    f << "\nVECTORS v float\n";
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                const vec &v = val(i, j, k).velocity();
                put<float>(f, v.x);
                put<float>(f, v.y);
                put<float>(f, v.z);
            }
    f << "\nSCALARS p float\nLOOKUP_TABLE default\n";
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                put<float>(f, gas().pressure(val(i, j, k)));
    f << "\nSCALARS T float\nLOOKUP_TABLE default\n";
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                put<float>(f, gas().temperature(val(i, j, k)));
    f << "\nSCALARS eps float\nLOOKUP_TABLE default\n";
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                put<float>(f, gas().specific_energy(val(i, j, k)));
    f << std::endl;
    f.close();
}

}

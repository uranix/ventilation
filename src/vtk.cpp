#include "../include/object.h"
#include <fstream>

namespace objects {

struct Appender {
    std::fstream &f;
    size_t written;
    Appender(std::fstream &f) : f(f), written(0) { }
    template<typename T>
    void putsize(size_t cnt) {
        put<int>(cnt * sizeof(T));
    }
    template<typename T>
    void put(T value) {
        f.write(reinterpret_cast<char *>(&value), sizeof(T));
        written += sizeof(T);
    }
};

struct Offsetter {
    size_t offset;
    Offsetter() : offset(0) { }
    template<typename T>
    size_t advance(size_t count) {
        size_t ret = offset;
        offset += sizeof(int) + count * sizeof(T);
        return ret;
    }
};

void object::save(const std::string &prefix, const int step) const {
    #include "GitVersion.h"
    std::string fn(prefix + id);
    fn += ".";
    fn += std::to_string(step);
    fn += ".vtr";

    std::fstream f(fn, std::ios::out | std::ios::binary);

    if (!f) {
        std::cerr << "Could not save file " << fn << std::endl;
        return;
    }

    Offsetter offs;
    Appender app(f);

    size_t ncells = nx * ny * nz;
    (void)ncells;

    f << "<!-- Written by: " << VERSION << " -->\n";
    f << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "  <RectilinearGrid WholeExtent = \"0 " << nx << " 0 " << ny << " 0 " << nz << "\">\n";
    f << "    <Piece Extent = \"0 " << nx << " 0 " << ny << " 0 " << nz << "\">\n";
    f << "      <PointData />\n";
    f << "      <CellData>\n";
    f << "        <DataArray type=\"Float32\" Name=\"v\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offs.advance<float>(3 * ncells) << "\"/>\n";
    for (int ic = 0; ic < nc; ic++)
        f << "        <DataArray type=\"Float32\" Name=\"density_" << ic << "\" format=\"appended\" offset=\"" << offs.advance<float>(ncells) << "\"/>\n";
    f << "        <DataArray type=\"Float32\" Name=\"density\" format=\"appended\" offset=\"" << offs.advance<float>(ncells) << "\"/>\n";

    f << "        <DataArray type=\"Float32\" Name=\"adiabatic_index\" format=\"appended\" offset=\"" << offs.advance<float>(ncells) << "\"/>\n";
    f << "        <DataArray type=\"Float32\" Name=\"pressure\" format=\"appended\" offset=\"" << offs.advance<float>(ncells) << "\"/>\n";
    f << "        <DataArray type=\"Float32\" Name=\"temperature\" format=\"appended\" offset=\"" << offs.advance<float>(ncells) << "\"/>\n";
    f << "        <DataArray type=\"Float32\" Name=\"specific_energy\" format=\"appended\" offset=\"" << offs.advance<float>(ncells) << "\"/>\n";
#if TURBULENCE
    f << "        <DataArray type=\"Float32\" Name=\"turb_k\" format=\"appended\" offset=\"" << offs.advance<float>(ncells) << "\"/>\n";
    f << "        <DataArray type=\"Float32\" Name=\"turb_eps\" format=\"appended\" offset=\"" << offs.advance<float>(ncells) << "\"/>\n";
    f << "        <DataArray type=\"Float32\" Name=\"turb_viscosity\" format=\"appended\" offset=\"" << offs.advance<float>(ncells) << "\"/>\n";
    f << "        <DataArray type=\"Float32\" Name=\"turb_production\" format=\"appended\" offset=\"" << offs.advance<float>(ncells) << "\"/>\n";
#endif
    f << "      </CellData>\n";
    f << "      <Coordinates>\n";
    f << "        <DataArray type=\"Float32\" Name=\"X_coordinates\" format=\"appended\" offset=\"" << offs.advance<float>(nx + 1) << "\"/>\n";
    f << "        <DataArray type=\"Float32\" Name=\"Y_coordinates\" format=\"appended\" offset=\"" << offs.advance<float>(ny + 1) << "\"/>\n";
    f << "        <DataArray type=\"Float32\" Name=\"Z_coordinates\" format=\"appended\" offset=\"" << offs.advance<float>(nz + 1) << "\"/>\n";
    f << "      </Coordinates>\n";
    f << "    </Piece>\n";
    f << "  </RectilinearGrid>\n";
    f << "  <AppendedData encoding=\"raw\">\n  _";

    app.putsize<float>(3 * ncells);
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                const vec &v = val(i, j, k).velocity();
                app.put<float>(v.x);
                app.put<float>(v.y);
                app.put<float>(v.z);
            }
    for (int ic = 0; ic < nc; ic++) {
        app.putsize<float>(ncells);
        for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                    app.put<float>(val(i, j, k).rho[ic]);
    }

    app.putsize<float>(ncells);
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                app.put<float>(val(i, j, k).density());

    app.putsize<float>(ncells);
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                app.put<float>(gas().gamma_ratio(val(i, j, k)));

    app.putsize<float>(ncells);
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                app.put<float>(gas().pressure(val(i, j, k)));

    app.putsize<float>(ncells);
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                app.put<float>(gas().temperature(val(i, j, k)));

    app.putsize<float>(ncells);
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                app.put<float>(gas().specific_energy(val(i, j, k)));

#if TURBULENCE
    app.putsize<float>(ncells);
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                app.put<float>(val(i, j, k).turb_energy());

    app.putsize<float>(ncells);
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                app.put<float>(val(i, j, k).turb_dissipation());

    app.putsize<float>(ncells);
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                app.put<float>(val(i, j, k).turb_viscosity());

    app.putsize<float>(ncells);
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                app.put<float>(val(i, j, k).Pk);
#endif

    /* XYZ coordinates */
    app.putsize<float>(nx + 1);
    for (int i = 0; i <= nx; i++)
        app.put<float>(ll.x + i * h.x);
    app.putsize<float>(ny + 1);
    for (int j = 0; j <= ny; j++)
        app.put<float>(ll.y + j * h.y);
    app.putsize<float>(nz + 1);
    for (int k = 0; k <= nz; k++)
        app.put<float>(ll.z + k * h.z);

    f << "\n  </AppendedData>\n";
    f << "</VTKFile>\n";

    if (app.written != offs.offset) {
        std::cerr << "Data size written does not match with the offsets!" << std::endl;
        abort();
    }

    f.close();
}

}

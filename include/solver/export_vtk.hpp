#pragma once

#include "core/params.hpp"
#include "core/types.hpp"
#include "physics/markers.hpp"
#include "solver/openlb_setup.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

using namespace olb;
using namespace olb::descriptors;

struct VtiArrayView {
  std::string name;
  std::string type;
  int components = 1;
  const void* data = nullptr;
  uint32_t bytes = 0;
};

class ImageDataVtiWriter {
public:
  ImageDataVtiWriter(int nx, int ny) : _nx(nx), _ny(ny) { }

  void writePointData(const std::string& path,
                      const std::vector<VtiArrayView>& arrays) const
  {
    std::ofstream vti(path, std::ios::binary);
    if (!vti.is_open()) return;

    vti << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "  <ImageData WholeExtent=\"0 " << (_nx - 1) << " 0 " << (_ny - 1)
        << " 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n"
        << "    <Piece Extent=\"0 " << (_nx - 1) << " 0 " << (_ny - 1) << " 0 0\">\n"
        << "      <PointData>\n";

    uint32_t offset = 0;
    for (const auto& a : arrays) {
      vti << "        <DataArray type=\"" << a.type
          << "\" Name=\"" << a.name << "\"";
      if (a.components > 1) {
        vti << " NumberOfComponents=\"" << a.components << "\"";
      }
      vti << " format=\"appended\" offset=\"" << offset << "\"/>\n";
      offset += 4 + a.bytes;
    }

    vti << "      </PointData>\n"
        << "    </Piece>\n"
        << "  </ImageData>\n"
        << "  <AppendedData encoding=\"raw\">\n_";

    for (const auto& a : arrays) {
      vti.write(reinterpret_cast<const char*>(&a.bytes), 4);
      vti.write(reinterpret_cast<const char*>(a.data), a.bytes);
    }

    vti << "\n  </AppendedData>\n</VTKFile>\n";
  }

private:
  int _nx;
  int _ny;
};

struct ExportSnapshot {
  int nx = 0;
  int ny = 0;
  int nPts = 0;
  std::vector<T> ux;
  std::vector<T> uy;
  std::vector<T> velocity;
  std::vector<T> uMag;
  std::vector<T> vorticity;
  std::vector<T> density;
  std::vector<T> forceMag;
  std::vector<int32_t> bodyMask;

  ExportSnapshot(int nx_, int ny_) : nx(nx_), ny(ny_), nPts(nx_ * ny_) {
    ux.resize(nPts);
    uy.resize(nPts);
    velocity.resize(3 * nPts);
    uMag.resize(nPts);
    vorticity.resize(nPts);
    density.resize(nPts);
    forceMag.resize(nPts);
    bodyMask.resize(nPts);
  }

  void sample(SuperLattice<T, DESCRIPTOR>& sLattice,
              const std::vector<int>& mask,
              bool includeDiagnostics)
  {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticeVelocity2D<T, DESCRIPTOR> velocityF(sLattice);
    SuperLatticeDensity2D<T, DESCRIPTOR> densityF(sLattice);
    SuperLatticeField2D<T, DESCRIPTOR, descriptors::FORCE> forceF(sLattice);

    for (int idx = 0; idx < nPts; ++idx) {
      const int j = idx / nx;
      const int i = idx % nx;
      int input[3] = {0, i, j};

      T u[2] = {0.0, 0.0};
      velocityF(u, input);
      ux[idx] = u[0];
      uy[idx] = u[1];

      if (mask[idx]) {
        velocity[3 * idx] = 0.0;
        velocity[3 * idx + 1] = 0.0;
        velocity[3 * idx + 2] = 0.0;
        uMag[idx] = 0.0;
      } else {
        velocity[3 * idx] = u[0];
        velocity[3 * idx + 1] = u[1];
        velocity[3 * idx + 2] = 0.0;
        uMag[idx] = std::sqrt(u[0] * u[0] + u[1] * u[1]);
      }

      if (includeDiagnostics) {
        T rho[1] = {1.0};
        T force[2] = {0.0, 0.0};
        densityF(rho, input);
        forceF(force, input);
        density[idx] = mask[idx] ? T(1) : rho[0];
        forceMag[idx] = std::sqrt(force[0] * force[0] + force[1] * force[1]);
        bodyMask[idx] = static_cast<int32_t>(mask[idx]);
      }
    }

    computeVorticity(mask);
  }

private:
  void computeVorticity(const std::vector<int>& mask)
  {
    std::fill(vorticity.begin(), vorticity.end(), T(0));
    for (int j = 1; j < ny - 1; ++j) {
      for (int i = 1; i < nx - 1; ++i) {
        const int idx = j * nx + i;
        if (mask[idx]) continue;
        vorticity[idx] = 0.5 * (uy[j * nx + (i + 1)] - uy[j * nx + (i - 1)])
                       - 0.5 * (ux[(j + 1) * nx + i] - ux[(j - 1) * nx + i]);
      }
    }
  }
};



inline std::vector<int> computeBodyMask(const LagrangianMarkers& mk, int nx, int ny)
{
  const int nPoly = mk.size();
  std::vector<int> mask(nx * ny, 0);
  if (nPoly < 3) return mask;
  T minPX = *std::min_element(mk.x.begin(), mk.x.end());
  T maxPX = *std::max_element(mk.x.begin(), mk.x.end());
  T minPY = *std::min_element(mk.y.begin(), mk.y.end());
  T maxPY = *std::max_element(mk.y.begin(), mk.y.end());
  int iLo = std::max(0,    static_cast<int>(std::floor(minPX)) - 1);
  int iHi = std::min(nx-1, static_cast<int>(std::ceil(maxPX))  + 1);
  int jLo = std::max(0,    static_cast<int>(std::floor(minPY)) - 1);
  int jHi = std::min(ny-1, static_cast<int>(std::ceil(maxPY))  + 1);
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j = jLo; j <= jHi; ++j) {
    for (int i = iLo; i <= iHi; ++i) {
      T px = static_cast<T>(i);
      T py = static_cast<T>(j);
      bool inside = false;
      for (int k = 0, l = nPoly - 1; k < nPoly; l = k++) {
        if (((mk.y[k] > py) != (mk.y[l] > py)) &&
            (px < (mk.x[l] - mk.x[k]) * (py - mk.y[k])
                  / (mk.y[l] - mk.y[k]) + mk.x[k]))
          inside = !inside;
      }
      if (inside) mask[j * nx + i] = 1;
    }
  }
  return mask;
}

inline void writeVelocityVTI(const std::string& velDir, int frame,
                      const ImageDataVtiWriter& vtiWriter,
                      const ExportSnapshot& fields)
{
  std::string fname = velDir + "eel3dof_velocity_" + std::to_string(frame) + ".vti";
  vtiWriter.writePointData(fname, {
    {"velocity", "Float64", 3, fields.velocity.data(),
     static_cast<uint32_t>(fields.velocity.size() * sizeof(T))},
    {"u_mag", "Float64", 1, fields.uMag.data(),
     static_cast<uint32_t>(fields.uMag.size() * sizeof(T))}
  });
}

inline void writeVorticityVTI(const std::string& vortDir, int frame,
                       const ImageDataVtiWriter& vtiWriter,
                       const ExportSnapshot& fields)
{
  std::string fname = vortDir + "eel3dof_vorticity_" + std::to_string(frame) + ".vti";
  vtiWriter.writePointData(fname, {
    {"vorticity", "Float64", 1, fields.vorticity.data(),
     static_cast<uint32_t>(fields.vorticity.size() * sizeof(T))}
  });
}

inline void writeDiagnosticsVTI(const std::string& diagDir, int frame,
                         const ImageDataVtiWriter& vtiWriter,
                         const ExportSnapshot& fields)
{
  std::string fname = diagDir + "eel3dof_diagnostics_" + std::to_string(frame) + ".vti";
  vtiWriter.writePointData(fname, {
    {"density", "Float64", 1, fields.density.data(),
     static_cast<uint32_t>(fields.density.size() * sizeof(T))},
    {"body_mask", "Int32", 1, fields.bodyMask.data(),
     static_cast<uint32_t>(fields.bodyMask.size() * sizeof(int32_t))},
    {"ibm_force_mag", "Float64", 1, fields.forceMag.data(),
     static_cast<uint32_t>(fields.forceMag.size() * sizeof(T))}
  });
}

inline void writeBodySnapshotCSV(const std::string& bodyOutDir, int frame,
                          const LagrangianMarkers& frameMarkers)
{
  std::string fname = bodyOutDir + "eel3dof_body_" + std::to_string(frame) + ".csv";
  std::ofstream csv(fname);
  if (!csv.is_open()) return;
  csv << "x,y\n";
  csv << std::setprecision(17);
  for (int m = 0; m < frameMarkers.size(); ++m)
    csv << frameMarkers.x[m] << "," << frameMarkers.y[m] << "\n";
}

inline void writeBodyVTP(const std::string& bodyOutDir, int frame,
                  const EelParams& p, const LagrangianMarkers& mk)
{
  const int nContour = mk.size();
  const int nSp = p.nSpine;
  const T dsSpine = p.eelScale * p.eelLength / (nSp - 1);
  const int nCap = std::max(static_cast<int>(std::ceil(M_PI * p.bodyRadius / dsSpine)), 8);
  const int nExpected = 2 * nSp + 2 * (nCap - 1);
  if (nContour != nExpected || nContour < 6) return;

  const int lowBase     = nSp + (nCap - 1);
  const int headCapBase = 2 * nSp + (nCap - 1);
  auto lowerMatch = [lowBase, nSp](int i) { return lowBase + (nSp - 1 - i); };

  const int nBodyTris = 2 * (nSp - 1);
  const int nTailTris = nCap - 1;
  const int nHeadTris = nCap - 1;
  const int nTris = nBodyTris + nTailTris + nHeadTris;

  std::string fname = bodyOutDir + "eel3dof_body_" + std::to_string(frame) + ".vtp";
  std::ofstream vtp(fname);
  if (!vtp.is_open()) return;
  vtp << std::setprecision(10);
  vtp << "<?xml version=\"1.0\"?>\n";
  vtp << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  vtp << "  <PolyData>\n";
  vtp << "    <Piece NumberOfPoints=\"" << nContour
      << "\" NumberOfVerts=\"0\" NumberOfLines=\"1\""
      << " NumberOfStrips=\"0\" NumberOfPolys=\"" << nTris << "\">\n";

  vtp << "      <Points>\n";
  vtp << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (int m = 0; m < nContour; ++m)
    vtp << "          " << mk.x[m] << " " << mk.y[m] << " 0\n";
  vtp << "        </DataArray>\n";
  vtp << "      </Points>\n";

  vtp << "      <Lines>\n";
  vtp << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  vtp << "         ";
  for (int m = 0; m < nContour; ++m) vtp << " " << m;
  vtp << " 0\n";
  vtp << "        </DataArray>\n";
  vtp << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  vtp << "          " << (nContour + 1) << "\n";
  vtp << "        </DataArray>\n";
  vtp << "      </Lines>\n";

  vtp << "      <Polys>\n";
  vtp << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for (int i = 0; i < nSp - 1; ++i) {
    int u0 = i, u1 = i + 1;
    int l0 = lowerMatch(i), l1 = lowerMatch(i + 1);
    vtp << "          " << u0 << " " << u1 << " " << l1 << "\n";
    vtp << "          " << u0 << " " << l1 << " " << l0 << "\n";
  }
  for (int k = 0; k < nCap - 1; ++k) {
    int a = nSp + k;
    int b = (k < nCap - 2) ? (nSp + k + 1) : lowBase;
    vtp << "          " << (nSp - 1) << " " << a << " " << b << "\n";
  }
  int headFan = lowBase + nSp - 1;
  for (int k = 0; k < nCap - 1; ++k) {
    int a = headCapBase + k;
    int b = (k < nCap - 2) ? (headCapBase + k + 1) : 0;
    vtp << "          " << headFan << " " << a << " " << b << "\n";
  }
  vtp << "        </DataArray>\n";
  vtp << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (int m = 0; m < nTris; ++m)
    vtp << "          " << 3 * (m + 1) << "\n";
  vtp << "        </DataArray>\n";
  vtp << "      </Polys>\n";

  vtp << "    </Piece>\n";
  vtp << "  </PolyData>\n";
  vtp << "</VTKFile>\n";
}

inline std::string writePVD(const std::string& runOutDir, const std::string& pvdName,
                     const std::string& relDir, const std::string& prefix,
                     const std::string& ext, int nFrames,
                     int exportInterval, T dtAnim)
{
  std::string pvdPath = runOutDir + pvdName;
  std::ofstream pvd(pvdPath);
  if (!pvd.is_open()) return pvdPath;
  pvd << "<?xml version=\"1.0\"?>\n";
  pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  pvd << "  <Collection>\n";
  for (int i = 0; i < nFrames; i += exportInterval) {
    T tPhys = (i + 1) * dtAnim;
    pvd << "    <DataSet timestep=\"" << tPhys
        << "\" file=\"" << relDir << prefix << i << ext << "\"/>\n";
  }
  pvd << "  </Collection>\n";
  pvd << "</VTKFile>\n";
  return pvdPath;
}

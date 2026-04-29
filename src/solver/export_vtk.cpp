// OpenLB's olb.h emits non-inline definitions in this build, so OpenLB-facing
// export implementations live inline in solver/export_vtk.hpp and are
// instantiated only by src/main.cpp.  This file intentionally remains empty
// to keep the target architecture stable without modifying ../openlb.

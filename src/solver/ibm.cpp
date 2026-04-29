// OpenLB-facing IBM code is intentionally implemented in solver/ibm.hpp.
// ibmStep is templated on the reset-force callable and the OpenLB lattice type;
// keeping it header-only avoids template link failures and avoids introducing
// another translation unit that includes olb.h in this OpenLB sample build.

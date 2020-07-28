#ifndef WRITE_TO_POSCAR_H
#define WRITE_TO_POSCAR_H
#include "interstitial_mesh.hpp"
#include "../avdv-factor-group/lattice.hpp"
#include "../avdv-factor-group/site.hpp"
#include "../avdv-factor-group/structure.hpp"
#include <fstream>
#include <string>
void write_to_poscar(const Structure& structure, std::ofstream& poscar_stream);


#endif

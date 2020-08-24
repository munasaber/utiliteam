#ifndef FACTOR_GROUP_HH
#define FACTOR_GROUP_HH

#include "./vectorfunctions.hpp"
#include "./symop.hpp"
#include "./site.hpp"
#include "./lattice.hpp"
#include "./structure.hpp"
#include "./symgroup.hpp"
#include "./point_group.hpp"

std::vector<Site> transform_basis(const SymOp& symop, const std::vector<Site>& basis);
bool basis_maps_onto_itself(const std::vector<Site>& original_basis, const std::vector<Site>& transformed_basis, const Lattice& lattice, double tol);

SymGroup<SymOp,BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f> generate_factor_group(const Structure& struc, double tol);

std::vector<Eigen::Vector3d> generate_translations(const Site& original_basis_site, const std::vector<Site>& transformed_basis);

#endif

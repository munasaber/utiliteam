#ifndef INTERSTITIAL_MESH_H
#define INTERSTITIAL_MESH_H

#include "../avdv-factor-group/lattice.hpp"
#include "../avdv-factor-group/vectorfunctions.hpp"
#include "../avdv-factor-group/site.hpp"
#include "../avdv-factor-group/symop.hpp"
#include "../avdv-factor-group/structure.hpp"
#include "../avdv-factor-group/symgroup.hpp"
#include <vector>
//Next tools we need:
//Geometric center of mass for cluster
//Find asymetric unit of a basis*
//Find all sites within a radius*

class Site;
class SymOp;

std::vector<Eigen::Vector3d> make_grid_points(int in_a, int in_b, int in_c, const Lattice& lattice);
//double distance(Coordinate& middlepoint, Coordinate& compared_coord);
//std::vector<Coordinate> find_sites_within_radius(Coordinate middlepoint, double my_radius, Structure& my_struc, Coordinate potential_anion_coord);
std::vector<Eigen::Vector3d> find_interstitials_within_radius(std::vector<Eigen::Vector3d>& interstitial_coordinates, const Eigen::Vector3d& sphere_origin, double radius, const Lattice& lattice);
std::vector<Eigen::Vector3d> keep_reasonable_interstitial_gridpoints(const std::vector<Eigen::Vector3d>& total_interstitial_coordinates, const std::vector<Eigen::Vector3d>& interstitial_coordinates_to_discard, double tol, Lattice& lattice);

std::vector<Eigen::Vector3d> make_orbit(const Eigen::Vector3d& coordinates, const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group, const Lattice& lattice);
std::vector<std::vector<Eigen::Vector3d>>
more_complex_bin_into_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  //const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group,
                                  const std::vector<SymOp>& symops,
		                  const Lattice& lattice,
                                  double tol);
std::vector<int>
label_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  //const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group,
                                  const std::vector<SymOp>& symops,
                                  const Lattice& lattice,
                                  double tol);
std::vector<std::vector<Eigen::Vector3d>>
bin_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  //const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group,
                                  const std::vector<SymOp>& symops,
                                  const Lattice& lattice,
                                  double tol);
std::vector<Eigen::Vector3d> make_asymmetric_unit(const std::vector<Eigen::Vector3d>& complete_structure_basis, 
                                  //const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group,
                                  const std::vector<SymOp>& symops,
				  const Lattice& lattice, 
			          double tol);
#endif 

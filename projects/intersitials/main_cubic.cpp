#include "interstitial_mesh.hpp"
#include "../avdv-factor-group/site.hpp"
#include "../avdv-factor-group/vectorfunctions.hpp"
#include "write_to_poscar.hpp"
#include "../avdv-factor-group/io.hpp"
#include "../avdv-factor-group/factor_group.hpp"
#include "../avdv-factor-group/point_group.hpp"
#include "../avdv-factor-group/lattice.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


	// auto point_group=generate_point_group(lattice, tol);
	//In the future consider doing list of symops instead of symgroup for interstitial_mesh stuff
int main()
{
	double tol=1e-5;
	Lattice lattice(Eigen::Vector3d(4,0,0), Eigen::Vector3d(0, 4, 0), Eigen::Vector3d(0, 0, 4));
	std::vector<Eigen::Vector3d> original_gridpoints=make_grid_points(10, 10, 10, lattice);
	auto point_group=generate_point_group(lattice,tol);
        std::vector<SymOp> symops=point_group.operations();
	auto  orbit_mesh=bin_by_symmetrical_equivalence(original_gridpoints, symops, lattice, tol);
	std::ofstream cubic_outputfile;
	int i=0;
	for (const auto& orbit: orbit_mesh)
	{
		std::string front("outputfiles/cubic_lattice_");
		std::string base(".vasp");
		std::vector<Site> sites_to_push_back;
		//cubic_outputfile.open();
			for (const auto& interstitial: orbit)
			{
				sites_to_push_back.emplace_back(Site("Li", (interstitial)));
			i++;
			}
		Structure cubic_with_sites(lattice, sites_to_push_back);
		cubic_outputfile.open(front+std::to_string(i)+base);
		write_to_poscar(cubic_with_sites, cubic_outputfile);

	}
	auto asym_unit=make_asymmetric_unit(original_gridpoints, symops, lattice, tol);
	std::cout<<asym_unit.size()<<std::endl;
}


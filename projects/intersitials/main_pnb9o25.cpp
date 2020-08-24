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


int main()
{
	double tol=1e-5;
	double radius=1;
	double orbit_tol=0.3;

	//PNB9O25
        Structure pnb9o25 = read_poscar("../avdv-factor-group/test_files/pnb9o25.vasp");
	std::vector<Eigen::Vector3d> original_gridpoints =make_grid_points(8, 16, 16, pnb9o25.get_lattice());

	std::vector<Eigen::Vector3d> coordinate_removal_list;
	std::vector<Eigen::Vector3d> original_atom_eigen_coordinates;
	for (int i=0; i<pnb9o25.get_sites().size(); i++)
	{
		original_atom_eigen_coordinates.push_back(pnb9o25.get_sites().at(i).get_eigen_coordinate());
	}
		
	for (const Eigen::Vector3d& original_atom_site : original_atom_eigen_coordinates)
	{
		std::vector<Eigen::Vector3d> interstitials_to_remove_per_radius=find_interstitials_within_radius(original_gridpoints, original_atom_site, radius, pnb9o25.get_lattice());	
		//do within radius for all atoms then compare within radius atoms to all atoms
		for (const Eigen::Vector3d& removal_atom: interstitials_to_remove_per_radius)
		{
			VectorPeriodicCompare_f test_coord(removal_atom, tol, pnb9o25.get_lattice());
			//VectorPeriodicCompare_f test_coord(removal_atom, tol);
			if(std::find_if(original_gridpoints.begin(), original_gridpoints.end(), test_coord)!=original_gridpoints.end())
			{
				coordinate_removal_list.push_back(removal_atom);
			
			}
		}
	}

	Lattice lattice=pnb9o25.get_lattice();
	std::vector<Eigen::Vector3d> remaining_interstitials=keep_reasonable_interstitial_gridpoints(original_gridpoints, coordinate_removal_list, tol, lattice);	
	
	//check what my outputs are for remaining_interstitials
	//for (const auto& interstital: remaining_interstitials)
	//{
	//	std::cout<< interstital<< std::endl;
	//}
	
	//get orbits so that you can print both full structure and  all the different orbit containers
	auto factor_group=generate_factor_group(pnb9o25, tol); 
        std::vector<SymOp> symops=factor_group.operations();
	std::ofstream orbitoutputfile;
	std::ofstream completeorbitoutputfile;
	std::vector<std::vector<Eigen::Vector3d>> orbit_container_list=bin_by_symmetrical_equivalence(remaining_interstitials, symops, lattice, orbit_tol);
	int i=0;
	std::vector<Site> complete_orbit_counter;
	for (const auto& orbit_container: orbit_container_list)
	{

		if (orbit_container.size()==4)
		{
			
        		std::string front("outputfiles/pnb9o25_");
			std::string base(".vasp");
			std::vector<Site> all_site_plus_indiv_orbit=pnb9o25.get_sites();
			for(const auto& interstitial_in_orbit_vector: orbit_container)
			{
				all_site_plus_indiv_orbit.emplace_back("Li", (interstitial_in_orbit_vector));
				complete_orbit_counter.emplace_back("Li", (interstitial_in_orbit_vector));
				i++;
			}
			orbitoutputfile.open(front+std::to_string(i)+base);
			Structure complete_structure_plus_indiv_orbit(lattice, all_site_plus_indiv_orbit);
			write_to_poscar(complete_structure_plus_indiv_orbit, orbitoutputfile);
			orbitoutputfile.close();
			
		}
	        if (orbit_container.size()>4)
		{
			std::cout<<"For some reason I have an orbit of size "<< orbit_container.size()<<std::endl;	
		}		

	}
	//get the complete set of all relevant orbits
	std::vector<Site> orbit_counter_plus_pnb9o25=pnb9o25.get_sites();
	for (const auto& orbit: complete_orbit_counter)
	{
		orbit_counter_plus_pnb9o25.emplace_back(orbit);
	}
	completeorbitoutputfile.open("outputfiles/pnb9o25_allorbitssize4.vasp");
	write_to_poscar(Structure(lattice, orbit_counter_plus_pnb9o25), completeorbitoutputfile);
	completeorbitoutputfile.close();

	std::cout<<orbit_container_list.size()<<std::endl;
	
	//std::cout<<asym_container.size();
        
//	std::ofstream asymoutputfile;
	//get 4 atom orbits from asym unit
	std::vector<Eigen::Vector3d> asym_container=make_asymmetric_unit(remaining_interstitials, symops, lattice, tol);
	std::cout<<asym_container.size();
//	int j=0;
//	for (const Eigen::Vector3d& asym_vector: asym_container)
//	{
//        	std::string front("outputfiles/pnb9o25_from_aym_");
//		std::string base(".vasp");
//	        std::vector<Site> four_sites_to_add_container;
//		four_sites_to_add_container.emplace_back("Li", Coordinate(asym_vector));
//		for (const SymOp& operation: factor_group.operations())
//		{
//			Eigen::Vector3d transformed_asym_vector=operation*asym_vector;
//			four_sites_to_add_container.emplace_back("Li", Coordinate(transformed_asym_vector));
//		}
//		asymoutputfile.open(front+std::to_string(j)+base);
//		Structure complete_structure_from_asym(lattice, four_sites_to_add_container);
//		write_to_poscar(complete_structure_from_asym, asymoutputfile);
//		j++;
//		
//	}
	
	//make new structure
	std::vector<Site> all_site_and_interstitials_container=pnb9o25.get_sites();
	for (const auto& eigen_coord: remaining_interstitials)
	{
		all_site_and_interstitials_container.emplace_back("Li", (eigen_coord));
	}
	Structure complete_structure(lattice, all_site_and_interstitials_container);
	
	std::ofstream completestructureoutputfile;
        std::string front("outputfiles/pnb9o25_completestructure");
	std::string base(".vasp");
	completestructureoutputfile.open(front+std::to_string(i)+base);
	//write_to_poscar(complete_structure, std::cout);
        write_to_poscar(complete_structure, completestructureoutputfile);
	return 0;
}


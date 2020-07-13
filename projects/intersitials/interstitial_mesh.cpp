#include "interstitial_mesh.hpp"
#include <algorithm>
#include <cmath>

// remove Coordinate and change to Eigen::Vector3d everywhere
std::vector<Eigen::Vector3d> make_grid_points(int in_a, int in_b, int in_c, const Lattice& lattice)
{
    // returns mxnxl grid of potential Li sites
    // convert to fractional...maybe
    std::vector<Eigen::Vector3d> grid;
    Eigen::VectorXd Divisions_in_a = Eigen::VectorXd::LinSpaced(in_a, 0, 1);
    Eigen::VectorXd Divisions_in_b = Eigen::VectorXd::LinSpaced(in_b, 0, 1);
    Eigen::VectorXd Divisions_in_c = Eigen::VectorXd::LinSpaced(in_c, 0, 1);
    // triple for loop to add in all the coordinate divisions
    for (int i = 0; i < in_a; i++)
    {
        for (int j = 0; j < in_b; j++)
        {
            for (int k = 0; k < in_c; k++)
            {
                Eigen::Vector3d temp_coord;
                temp_coord << Divisions_in_a(i), Divisions_in_b(j), Divisions_in_c(k);
                Eigen::Vector3d temp_coord_cart = convert_to_cartesian(lattice, temp_coord);
                grid.emplace_back(temp_coord_cart);
            }
        }
    }
    // make sure this comes out as cartesian
    return grid;
}

// returns the coordinates of interstitial sites within the radius of a certain atom. The intent is to take in intersitial/atom coords in
// cartesian  coords, but a radius in cartesian coords, and then return the interstitial sites within the radius in cartesian coords. Does
// it do this? Probably?
std::vector<Eigen::Vector3d> find_interstitials_within_radius(std::vector<Eigen::Vector3d>& interstitial_coordinates,
                                                              const Eigen::Vector3d& sphere_origin,
                                                              double radius,
                                                              const Lattice& lattice)
{

    // subtract sphere origin and then do distance algorithim (ie use peridic compare)
    // use radius as the tolerarnace????
    std::vector<Eigen::Vector3d> interstitials_within_radius;
    for (const auto& interstitial_coord : interstitial_coordinates)
    {
        // can probably do this without normalization
        Eigen::Vector3d normalized_interstitial_coord = interstitial_coord - sphere_origin;
        VectorPeriodicCompare_f compare_interstitials_relative_to_single_site(Eigen::Vector3d(0, 0, 0), radius, lattice);
        if (compare_interstitials_relative_to_single_site(normalized_interstitial_coord) == true)
        {
            interstitials_within_radius.emplace_back(interstitial_coord);
        }
    }
    return interstitials_within_radius;
}

// Discard certain interstital grid points within certain radius
std::vector<Eigen::Vector3d> keep_reasonable_interstitial_gridpoints(const std::vector<Eigen::Vector3d>& total_interstitial_coordinates, const std::vector<Eigen::Vector3d>& interstitial_coordinates_to_discard, double tol, Lattice& lattice)
{
	std::vector<Eigen::Vector3d> remaining_interstitials;
	for (const Eigen::Vector3d& interstitial_coord : total_interstitial_coordinates)
	{
		VectorPeriodicCompare_f compare_total_to_discarding_interstitials(interstitial_coord, tol, lattice);
		if (std::find_if(interstitial_coordinates_to_discard.begin(), interstitial_coordinates_to_discard.end(),
 compare_total_to_discarding_interstitials)==interstitial_coordinates_to_discard.end())
		{
			remaining_interstitials.emplace_back(interstitial_coord);
		}
	}
	return remaining_interstitials;
}

// TODO:

std::vector<Eigen::Vector3d> make_orbit(const Eigen::Vector3d& coordinate,
                                        const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group,
                                        const Lattice& lattice)
{
	std::vector<Eigen::Vector3d> orbit;
        for (const SymOp& symop : factor_group.operations())
        {
            Eigen::Vector3d transformedcoord = symop * coordinate;
    	    orbit.emplace_back(transformedcoord);	    
        }
	return orbit;
}
// Indexes vectors of the different asymmetric orbits of the system. This cam further be broken into the asymmetric atoms in the system
// For reference these use the interstitial coordinates after the old ones are taken out
// The definition for orbit is vague and confusing and I do not like it

std::vector<std::vector<Eigen::Vector3d>>
more_complex_bin_into_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group,
                                  const Lattice& lattice,
                                  double tol)
{
    std::vector<std::vector<Eigen::Vector3d>> orbit_container;
    for (const auto& interstitial_coordinate : coordinates)
    {
        for (const SymOp& symop : factor_group.operations())
        {
	    bool already_found=false;
	    Eigen::Vector3d transformedcoord = symop * interstitial_coordinate;
            for (std::vector<Eigen::Vector3d>& orbit : orbit_container)
            {
                VectorPeriodicCompare_f test_coord(transformedcoord, tol, lattice);
                if (find_if(orbit.begin(), orbit.end(), test_coord) != orbit.end())
                {
                     already_found=true;
		     //consider doing vector periodic compare for interstitial is something I should keep
		     VectorPeriodicCompare_f make_sure_interstitial_is_unique(interstitial_coordinate, tol, lattice);
		     if (find_if(orbit.begin(), orbit.end(), make_sure_interstitial_is_unique)==orbit.end())
		     	orbit.emplace_back(interstitial_coordinate);
                }
            }
            // if the site is not in any of the existing orbits then we need to create
            // a new orbit
            if (!already_found)
            {
                std::vector<Eigen::Vector3d> new_orbit;
                new_orbit.emplace_back(interstitial_coordinate);
		orbit_container.emplace_back(new_orbit);
            }
        }
    }
    return orbit_container;
}

std::vector<int>
label_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group,
                                  const Lattice& lattice,
                                  double tol)
{
    std::vector<int> coordinate_tags(coordinates.size(),-1);
    for(int ix = 0; ix<coordinate_tags.size(); ++ix)
    {
        if(coordinate_tags[ix]!=-1)
        {
            continue;
        }

        const Eigen::Vector3d& coord=coordinates[ix];
        std::vector<Eigen::Vector3d> coord_orbit=make_orbit(coord,factor_group,lattice);

        for(int ixx=ix; ixx<coordinate_tags.size(); ++ixx)
        {
            if(coordinate_tags[ixx]!=-1)
            {
                continue;
            }

            VectorPeriodicCompare_f equals_ixx_coord(coordinates[ixx],tol,lattice);
            if(std::find_if(coord_orbit.begin(),coord_orbit.end(),equals_ixx_coord)!=coord_orbit.end())
            {
                coordinate_tags[ixx]=ix;
            }
        }
    }

    for(auto l : coordinate_tags)
    {
        assert(l!=-1);
    }

    return coordinate_tags;
}


std::vector<std::vector<Eigen::Vector3d>>
bin_by_symmetrical_equivalence(const std::vector<Eigen::Vector3d>& coordinates,
                                  const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group,
                                  const Lattice& lattice,
                                  double tol)
{
	std::vector<int> coordinate_tags= label_by_symmetrical_equivalence(coordinates, factor_group, lattice, tol);
	//goes through length of coordinate list and connects coordinate to orbit	
	int num_orbits=*std::max_element(coordinate_tags.begin(), coordinate_tags.end())+1;
	std::vector<std::vector<Eigen::Vector3d>> orbit_container;
	orbit_container.resize(num_orbits);
       	for (int i=0; i<coordinates.size(); i++)
	{
		Eigen::Vector3d temp_coord=coordinates[i];
		int label=coordinate_tags[i];
		VectorCompare_f compare_temp_coord_to_whats_in_index_of_orbit_container(temp_coord, tol);
		if (find_if(orbit_container[label].begin(), orbit_container[label].end(), compare_temp_coord_to_whats_in_index_of_orbit_container)==orbit_container[label].end())
		{
		orbit_container[label].emplace_back(temp_coord);
		}
	}
	return orbit_container;	
}
/// Tracks the assymetric intersttial atoms
//std::vector<Eigen::Vector3d> make_asymmetric_unit(const std::vector<Eigen::Vector3d>& complete_structure_basis,
//                                                  const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group,
//                                                  const Lattice& lattice,
//                                                  double tol)
//{
//    std::vector<Eigen::Vector3d> asymmetric_unit;
//    // for each site on the basis
//    // apply allsymmetry operations
//    for (const auto& basis : complete_structure_basis)
//    {
//        for (const auto& Symmetry_operation : factor_group.operations())
//        {
//            Eigen::Vector3d transformedcoord = Symmetry_operation * basis;
//            // TODO: Which comparator should you use?
//            VectorPeriodicCompare_f test_coord(transformedcoord, tol, lattice);
//
//            if (find_if(asymmetric_unit.begin(), asymmetric_unit.end(), test_coord) == asymmetric_unit.end())
//            {
//                asymmetric_unit.push_back(transformedcoord);
//            }
//        }
//    }
//    return asymmetric_unit;
//}
std::vector<Eigen::Vector3d> make_asymmetric_unit(const std::vector<Eigen::Vector3d>& complete_structure_basis,
                                                  const SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f>& factor_group,
                                                  const Lattice& lattice,
                                                  double tol)
{
	std::vector<Eigen::Vector3d> asymmetric_unit_collated;
	std::vector<std::vector<Eigen::Vector3d>> total_orbits=bin_by_symmetrical_equivalence(complete_structure_basis, factor_group, lattice, tol);
	for (const auto& orbit: total_orbits)
	{
		asymmetric_unit_collated.push_back(orbit[0]);
	}
	return asymmetric_unit_collated;
}

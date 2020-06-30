#include "../avdv-factor-group/coordinate.hpp"
#include "../avdv-factor-group/io.hpp"
#include "../avdv-factor-group/lattice.hpp"
#include "../avdv-factor-group/structure.hpp"
#include "../avdv-factor-group/tests.hpp"
#include "interstitial_mesh.hpp"

bool does_find_sites_within_radius_work_two_sites()
{
    double radius = 0.1;
    Eigen::Vector3d within_radius(Eigen::Vector3d(0.55, 0.55, 0.55));
    Eigen::Vector3d outside_radius(Eigen::Vector3d(0.9, 0.9, 0.9));
    std::vector<Eigen::Vector3d> all_interstitial_sites{within_radius, outside_radius};
    // I don't like this way of doing lattice at all.
    return find_interstitials_within_radius(all_interstitial_sites, Eigen::Vector3d(0.5, 0.5, 0.5), radius,
                                            Lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1)))
               .size() == 1;
}

bool does_find_sites_within_radius_work_for_exact_coordinates_general()
{
    double radius = 0.1;
    Eigen::Vector3d within_radius(Eigen::Vector3d(0.55, 0.55, 0.55));
    Eigen::Vector3d outside_radius(Eigen::Vector3d(0.9, 0.9, 0.9));
    std::vector<Eigen::Vector3d> all_interstitial_sites{within_radius, outside_radius};
    Eigen::Vector3d sphere_origin(Eigen::Vector3d(0.5, 0.5, 0.5));
    return find_interstitials_within_radius(all_interstitial_sites, sphere_origin, radius,
                                            Lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1)))
        .at(0)
        .isApprox(within_radius);
}

bool does_find_sites_within_radius_work_for_pnb9o25()
{
    double radius=1.95; //A
    Structure pnb9o25 = read_poscar("../avdv-factor-group/test_files/pnb9o25.vasp");
    // need to decide whether to do this based on cartesian coordinates or not? I'm assuming yes?
    Eigen::Vector3d sphere_origin(Eigen::Vector3d(0, 0, 0));
    Eigen::Vector3d within_radius_1_in_fractional(Eigen::Vector3d(-0.11500,  0.16110,  0.06890));
    Eigen::Vector3d within_radius_1=convert_to_cartesian(pnb9o25.get_lattice(), within_radius_1_in_fractional);
    Eigen::Vector3d within_radius_2_in_fractional(Eigen::Vector3d(0.67090,  0.21330,  0.44490));
    Eigen::Vector3d within_radius_2=convert_to_cartesian(pnb9o25.get_lattice(), within_radius_2_in_fractional);
    std::vector<Eigen::Vector3d> all_interstitial_sites{within_radius_1, within_radius_2};
    ////change this!!!!!

    return find_interstitials_within_radius(all_interstitial_sites, sphere_origin, radius, pnb9o25.get_lattice()).size() ==
           1;
}



bool does_find_sites_within_radius_work_for_pnb9o25_exact_coordinates()
{
    Structure pnb9o25 = read_poscar("../avdv-factor-group/test_files/pnb9o25.vasp");


    double radius=1.95; //A
    Eigen::Vector3d sphere_origin(Eigen::Vector3d(0, 0, 0));
    Eigen::Vector3d within_radius_1_in_fractional(Eigen::Vector3d(-0.11500,  0.16110,  0.06890));
    Eigen::Vector3d within_radius_1=convert_to_cartesian(pnb9o25.get_lattice(), within_radius_1_in_fractional);
    Eigen::Vector3d within_radius_2_in_fractional(Eigen::Vector3d(-0.5,  0.0,  0.0));
    Eigen::Vector3d within_radius_2=convert_to_cartesian(pnb9o25.get_lattice(), within_radius_2_in_fractional);
    std::vector<Eigen::Vector3d> all_interstitial_sites{within_radius_1, within_radius_2};

    std::vector<Eigen::Vector3d> coords_within_radius =
        find_interstitials_within_radius(all_interstitial_sites, sphere_origin, radius, pnb9o25.get_lattice());

    for (const Eigen::Vector3d expected_in_radius : {within_radius_1,within_radius_1})
    {
        VectorCompare_f test_coord(expected_in_radius,1e-5);

        if (find_if(coords_within_radius.begin(), coords_within_radius.end(), test_coord) == coords_within_radius.end())
        {
            return false;
        }
    }
    return true;
}

bool does_find_sites_within_radius_work_for_pnb9o25_all_corners_of_octahedra()
{

	Structure pnb9o25=read_poscar("../avdv-factor-group/test_files/pnb9o25.vasp");
	double radius=1.95;
	Eigen::Vector3d sphere_origin(Eigen::Vector3d(0,0,0));

	//this gets all the oxygen coordinates around the 0,0,0 octahedra and converts them to cartesian
	Eigen::Vector3d within_radius_1_in_fractional(Eigen::Vector3d(0.11500,  0.83890,  0.93110));
	Eigen::Vector3d within_radius_1=convert_to_cartesian(pnb9o25.get_lattice(), within_radius_1_in_fractional);
	Eigen::Vector3d within_radius_2_in_fractional(Eigen::Vector3d(-0.11500,  0.16110,  0.06890));	
	Eigen::Vector3d within_radius_2=convert_to_cartesian(pnb9o25.get_lattice(), within_radius_2_in_fractional);
	Eigen::Vector3d within_radius_3_in_fractional(Eigen::Vector3d(-0.50000,  0.00000,  0.00000));
	Eigen::Vector3d within_radius_3=convert_to_cartesian(pnb9o25.get_lattice(), within_radius_3_in_fractional);
	Eigen::Vector3d within_radius_4_in_fractional(Eigen::Vector3d(0.50000,  0.00000,  0.0000));	
	Eigen::Vector3d within_radius_4=convert_to_cartesian(pnb9o25.get_lattice(), within_radius_4_in_fractional);
	Eigen::Vector3d within_radius_5_in_fractional(Eigen::Vector3d(-0.04610, -0.06890,  0.16110));
	Eigen::Vector3d within_radius_5=convert_to_cartesian(pnb9o25.get_lattice(), within_radius_5_in_fractional);
	Eigen::Vector3d within_radius_6_in_fractional(Eigen::Vector3d(0.04610,  0.06890, -0.16110));
	Eigen::Vector3d within_radius_6=convert_to_cartesian(pnb9o25.get_lattice(), within_radius_6_in_fractional);
	std::vector<Eigen::Vector3d> all_within_radius{within_radius_1, within_radius_2, within_radius_3, within_radius_4, within_radius_5, within_radius_6};
	//this gets the coordinate of a far phosphorus atom
	Eigen::Vector3d outside_radius_in_fractional(Eigen::Vector3d(0.25000,  0.50000,  0.50000));
	Eigen::Vector3d outside_radius=convert_to_cartesian(pnb9o25.get_lattice(), outside_radius_in_fractional);
	std::vector<Eigen::Vector3d> all_atoms{within_radius_1, within_radius_2, within_radius_3, within_radius_4, within_radius_5, within_radius_6, outside_radius};
	std::vector<Eigen::Vector3d> coords_within_radius=find_interstitials_within_radius(all_atoms, sphere_origin, radius, pnb9o25.get_lattice());

	//check if all coordinates are correct, ie its not also inclusing the phosphorous as within
	for (const Eigen::Vector3d expected_in_radius :  all_within_radius)
	{
		VectorCompare_f test_coord(expected_in_radius, 1e-5);	
		
		if (find_if(coords_within_radius.begin(), coords_within_radius.end(), test_coord)==coords_within_radius.end())
		{
			return false;
		}
	}
	return true;
}




bool does_keep_reasonable_interstitial_gridpoints_work(double tol)
{
	Lattice unit_lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
	Eigen::Vector3d within_radius_1(Eigen::Vector3d(0.55, 0.55, 0.55));
	Eigen::Vector3d outside_radius_1(Eigen::Vector3d(0.9, 0.9, 0.9));
	Eigen::Vector3d within_radius_2(Eigen::Vector3d(0.74, 0.51, 0.51));
	Eigen::Vector3d outside_radius_2(Eigen::Vector3d(-0.24, 0.5, -0.5));
	std::vector<Eigen::Vector3d> within_radius{within_radius_1, within_radius_2};
	std::vector<Eigen::Vector3d> outside_radius{outside_radius_1, outside_radius_2};
	std::vector<Eigen::Vector3d> total_interstitials{within_radius_1, within_radius_2, outside_radius_1, outside_radius_2};
	return keep_reasonable_interstitial_gridpoints(total_interstitials, within_radius, tol, unit_lattice).size()==2;
	//return keep_reasonable_interstitial_gridpoints(total_interstitials, within_radius, tol, Lattice(Eigen::Vector3d(1,0,0), Eigen::Vector3d(0,1,0), Eigen::Vector3d(0,0,1))).size()==2; 
}



bool does_keep_reasonable_interstitial_gridpoints_work_for_exact_coordinates(double tol)
{
	
	Lattice unit_lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
	Eigen::Vector3d within_radius_1(Eigen::Vector3d(0.55, 0.55, 0.55));
	Eigen::Vector3d outside_radius_1(Eigen::Vector3d(0.9, 0.9, 0.9));
	Eigen::Vector3d within_radius_2(Eigen::Vector3d(0.74, 0.51, 0.51));
	Eigen::Vector3d outside_radius_2(Eigen::Vector3d(-0.24, 0.5, -0.5));
	std::vector<Eigen::Vector3d> within_radius{within_radius_1, within_radius_2};
	std::vector<Eigen::Vector3d> outside_radius{outside_radius_1, outside_radius_2};
	std::vector<Eigen::Vector3d> total_interstitials{within_radius_1, within_radius_2, outside_radius_1, outside_radius_2};
	std::vector<Eigen::Vector3d> remaining_gridpoints= keep_reasonable_interstitial_gridpoints(total_interstitials, within_radius, tol, unit_lattice);
 	for (const Eigen::Vector3d& expected_radius_coord : outside_radius)
	{
		VectorCompare_f test_coord(expected_radius_coord, tol);
		if (find_if(remaining_gridpoints.begin(), remaining_gridpoints.end(), test_coord)==remaining_gridpoints.end())
		{
			return false;
		}

	}
	return true;
}

bool does_make_orbit_work()
{
	return 0;
};


bool does_bin_into_symmetrically_equivalent_work_unit_lattice(double tol)
{
    Lattice unit_lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));

    Eigen::Vector3d base_coordinate(Eigen::Vector3d(0.51, 0.51, 0.75));
    Eigen::Vector3d symmetrically_equivalent(Eigen::Vector3d(0.51, 0.51, -0.75));
    Eigen::Vector3d symmetrically_ineqivalent(Eigen::Vector3d(0.9, 0.9, 0.9));

    SymOp mirror(Eigen::Matrix3d(1, 0, 0, 0, 1, 0, 0, 0, -1));
    BinarySymOpPeriodicCompare_f comparison(unit_lattice, tol);
    BinarySymOpPeriodicMultiplier_f mult_op(unit_lattice, tol);
    SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f> factor_group({mirror}, comparison, mult_op);
    std::vector<Eigen::Vector3d> all_interstitial_coordinates{base_coordinate, symmetrically_ineqivalent, symmetrically_equivalent};

    //TODO: Rename to bin into orbits, or bin into equivalents, or something.
    //so...why doesn't this work?
    return bin_into_symmetrically_equivalent(all_interstitial_coordinates, factor_group, unit_lattice, tol).size() == 2;
}

bool does_label_by_symmetrical_equivalence_work_unit_lattice(double tol)
{
    Lattice unit_lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
    Eigen::Vector3d base_coordinate(0.5, 0.5, 0.75);
    Eigen::Vector3d symmetrically_equivalent(Eigen::Vector3d(0.5001, 0.5001, -0.75));
    Eigen::Vector3d symmetrically_inequivalent(Eigen::Vector3d(0.75, 0.75, 1));
    SymOp mirror(Eigen::Matrix3d(1, 0, 0, 0, 1, 0, 0, 0, -1));
    BinarySymOpPeriodicCompare_f compare_coordinates(unit_lattice, tol);
    BinarySymOpPeriodicMultiplier_f mult_op(unit_lattice, tol);
    SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f> factor_group({mirror}, compare_coordinates, mult_op);   
    std::vector<Eigen::Vector3d> all_interstitial_coordinates{base_coordinate, symmetrically_inequivalent, symmetrically_equivalent};
    return label_by_symmetrical_equivalence(all_interstitial_coordinates, factor_group, unit_lattice, tol).size()==2;

}
bool does_make_orbits_work_for_containing_coordinates() { return 0; }

bool does_make_asymmetric_unit_work_for_pnb9o25(double tol)
{
    Structure pnb9o25 = read_poscar("../avdv-factor-group/test_files/pnb9o25.vasp");
    Eigen::Vector3d base_coordinate(0.25000,  0.50000,  0.50000);
    Eigen::Vector3d niobium_coordinate(
}

bool does_make_grid_points_work()
{
    int a1 = 5;
    int a2 = 3;
    int a3 = 4;
    if (make_grid_points(a1, a2, a3, Lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1))).size() == 60)
    {
        return true;
    }
    return false;
}

int main()
{
    double tol = 1e-6;
    double radius = 0.1;
    EXPECT_TRUE(does_find_sites_within_radius_work_two_sites(), "find interstitials sites within radius (should be 1)");
    EXPECT_TRUE(does_find_sites_within_radius_work_for_exact_coordinates_general(),
                "I get the correct coordinate that is within the radius");
    EXPECT_TRUE(does_find_sites_within_radius_work_for_pnb9o25(),
                "find interstitials sites within radius of site in pnb9o25 lattice (should be 2)");
    EXPECT_TRUE(does_find_sites_within_radius_work_for_pnb9o25_exact_coordinates(), "should have both correct coordinates in there");
    EXPECT_TRUE(does_find_sites_within_radius_work_for_pnb9o25_all_corners_of_octahedra(), "all octahedral corners (6) should be within the radius of the Niobium but P shouldn't");
    EXPECT_TRUE(does_keep_reasonable_interstitial_gridpoints_work(tol), "does keep reasonable interstitial gridpoints get the correct vector size (2)");
    EXPECT_TRUE(does_keep_reasonable_interstitial_gridpoints_work_for_exact_coordinates(tol), "does keep reasonable interstitial functions get the right coordinates");
    EXPECT_TRUE(does_bin_into_symmetrically_equivalent_work_unit_lattice(tol), "does bin into symmetrically equivalent orbits work for unit lattice");
    EXPECT_TRUE(does_label_by_symmetrical_equivalence_work_unit_lattice(tol), "does label by symmetricla equivalence work for a unit lattice");
    EXPECT_TRUE(does_make_grid_points_work(), "Check that I can appropriately make grid points");
}

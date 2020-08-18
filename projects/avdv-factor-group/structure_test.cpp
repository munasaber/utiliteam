#include "structure.hpp"
#include "tests.hpp"
#include <iostream>

int main()
{

    std::cout<<"---- Runnning Structure Tests ----"<<std::endl;
    std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
    //test constructor
std::vector<Site> test_sites;
Lattice my_lattice=Lattice(Eigen::Vector3d(1, 1, 0), Eigen::Vector3d(1, 0, 1), Eigen::Vector3d(0, 1, 1));
test_sites.push_back(Site("Na", Eigen::Vector3d(0.1, 0.1, 0.1)));
Structure my_structure=Structure(my_lattice,test_sites);

//test get lattice
EXPECT_TRUE(my_structure.get_lattice().row_vector_matrix()==my_lattice.row_vector_matrix(), "test Structure's get_lattice");

//test site compare

//test get sites
SiteCompare_f test_compare_sites(test_sites.at(0), 0.00001);
EXPECT_TRUE(test_compare_sites(my_structure.get_sites().at(0)), "test Structure's get_sites");
return 0;
}

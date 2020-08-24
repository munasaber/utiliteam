#include "symop.hpp"
#include "site.hpp"
#include "tests.hpp"
#include <iostream>

int main()
{
    std::cout<<"---- Running Site Tests ----"<<std::endl;
    std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
    //test constructor 
Site my_site= Site("Na", Eigen::Vector3d(0.1, 0.1, 0.1));

//test get atom
std::string test_atom="Na";
EXPECT_TRUE(my_site.get_atom()==test_atom, "check test_atom");

//test get_eigen_coordinate
Eigen::Vector3d test_eigen_coordinate;
test_eigen_coordinate<<0.1,0.1,0.1;
EXPECT_TRUE(my_site.get_eigen_coordinate()==test_eigen_coordinate, "check test_eigen_coordinate");

//test get_coordinate
EXPECT_TRUE(my_site.get_eigen_coordinate()==test_eigen_coordinate, "check test_eigen_coordinate");

//test site compare
SiteCompare_f my_site_compare(my_site, 0.0001);
Site test_site= Site("Na", (Eigen::Vector3d(0.1, 0.1, 0.1)));
EXPECT_TRUE(my_site_compare(test_site), "test site compare");
Site test_wrong_site= Site("Na", (Eigen::Vector3d(0.5, 0.1, 0.1)));
EXPECT_TRUE((!my_site_compare(test_wrong_site)), "test site compare but with the wrong function");

//site periodic compare
Lattice my_lattice=Lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 2, 0), Eigen::Vector3d(0, 1, 3));
SitePeriodicCompare_f my_site_periodic_compare(my_site, 0.0001, my_lattice);
Site test_site_periodic= Site("Na", (Eigen::Vector3d(1.1, 0.1, 0.1)));
EXPECT_TRUE(my_site_periodic_compare(test_site_periodic), "test site periodic compare");

//test periosic compare near boundaries
Lattice identity_lattice=Lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
Site test_site_close_to_1= Site("Na", (Eigen::Vector3d(0.9999999,0.9999999, 0.9999999)));
Site test_site_close_to_0= Site("Na", (Eigen::Vector3d(0.0000001, 0.0000001, 0.0000001)));
SitePeriodicCompare_f periodic_compare_close_to_zero(test_site_close_to_0, 0.0001, identity_lattice);
EXPECT_TRUE(periodic_compare_close_to_zero(test_site_close_to_1), "test site periodic compare close to boundaries");

//site periodic compare with negaitve value
Eigen::Matrix3d test_row_vector_lattice_neg;
Lattice lattice_neg=Lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 2, 0), Eigen::Vector3d(0, 1, 3));
SitePeriodicCompare_f my_site_periodic_compare_neg(my_site, 0.0001, my_lattice);
Site test_site_periodic_neg= Site("Na", (Eigen::Vector3d(-1.9, 0.1, 0.1)));
EXPECT_TRUE(my_site_periodic_compare_neg(test_site_periodic_neg), "test site periodic compare with negative value in x direction");

//operator*

Eigen::Matrix3d inversion_matrix;
inversion_matrix<<1,0,0,0,1,0,0,0,-1;
SymOp inversion =SymOp(inversion_matrix, Eigen::Vector3d(0,0,0));
Eigen::Vector3d unit_basis;
unit_basis<<1,1,1;
Eigen::Vector3d correct_basis;
correct_basis<<1,1,-1;
Lattice operator_lattice=Lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
SitePeriodicCompare_f site_compare_for_operator(Site("Na", correct_basis), 0.0001, operator_lattice);
Site comparing_site=Site("Na", unit_basis);
auto multiplication=inversion*comparing_site;
EXPECT_TRUE(site_compare_for_operator(multiplication), "testing operator function with default translation");


//DONE
return 0;
}

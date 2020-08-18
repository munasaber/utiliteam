#include "tests.hpp"
#include "lattice.hpp"
#include <iostream>

int main()
{

    std::cout<<"---- Running Lattice Tests ----"<<std::endl;
    std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
Eigen::Matrix3d test_row_vector_lattice;
test_row_vector_lattice<<1, 0, 0, 0, 2, 0, 0, 1, 3;
Eigen::Matrix3d test_col_vector_lattice=test_row_vector_lattice.transpose();
Lattice my_lattice=Lattice(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 2, 0), Eigen::Vector3d(0, 1, 3));
EXPECT_TRUE(my_lattice.row_vector_matrix()==test_row_vector_lattice, "row lattice test");
EXPECT_TRUE(my_lattice.col_vector_matrix()==test_col_vector_lattice, "col lattice test");
Eigen::Vector3d test_lattice_vector;
test_lattice_vector<<1, 0, 0;
EXPECT_TRUE(my_lattice.lattice_vector(0)==test_lattice_vector, "1st row vector test");

//test convert to fractional
Lattice doubled_lattice(Eigen::Vector3d(2, 0, 0), Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1));
Eigen::Vector3d fractional_vector;
fractional_vector<<1,0,0;
Eigen::Vector3d cartesian_vector;
cartesian_vector<<2,0,0;
EXPECT_TRUE(convert_to_fractional(doubled_lattice, cartesian_vector).isApprox(fractional_vector, 0.0001), "test conversion from cartesian to fractional");


//test bring within
Lattice unit_lattice(Eigen::Vector3d(1,0,0), Eigen::Vector3d(0,1,0),Eigen::Vector3d(0,0,1));
Eigen::Vector3d original_translation;
original_translation<<1.5, 0, 0;
Eigen::Vector3d brought_within_translation;
brought_within_translation<<0.5, 0, 0;
EXPECT_TRUE(bring_within(unit_lattice, 0.0001, original_translation).isApprox(brought_within_translation, 0.0001), "test the bring within function");


//test bring within near edge
Eigen::Vector3d original_translation_near_edge;
original_translation_near_edge<<-.001, 0, 0;
Eigen::Vector3d brought_within_translation_near_edge;
brought_within_translation_near_edge<<.999, 0, 0;
EXPECT_TRUE(bring_within(unit_lattice, 0.0001, original_translation_near_edge).isApprox(brought_within_translation_near_edge, 0.0001), "test the bring within function near edge");

//test bring within near corner
Eigen::Vector3d original_translation_near_corner;
original_translation_near_corner<<-.001, .999, 0;
Eigen::Vector3d brought_within_translation_near_corner;
brought_within_translation_near_corner<<.999, 0.999, 0;
EXPECT_TRUE(bring_within(unit_lattice, 0.0001, original_translation_near_corner).isApprox(brought_within_translation_near_corner, 0.0001), "test the bring within function near corner");
return 0;
}

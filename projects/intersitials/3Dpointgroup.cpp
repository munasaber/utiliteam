#include "xtal_classes.hpp"
// Initial Class set up
// Calculate the grid from which we will calculate potential Lprime values. Need 3 by 3 integers

std::vector<Eigen::Vector3d> calculate_gridpoints(Eigen::Matrix3d lattice,
                                                  int radius) // Might take in a structure instead and get a lattice from it
{
    // Eigen::Matrix3d lattice=my_lattice.col_vector_matrix();
    std::vector<Eigen::Vector3d> gridpoints; // left in int since I put radius in int
    Eigen::Vector3d pn;
    for (int m = -radius; m < (radius + 1); m++)
    {
        for (int n = -radius; n < (radius + 1); n++)
        {
            for (int j = -radius; j < (radius + 1); j++)
            {
                pn = n * lattice.col(0) + m * lattice.col(1) + j * lattice.col(2);
                gridpoints.push_back(pn);
            }
        }
    }

    return gridpoints;
}
// change complete

// calculate a list of potentials Lprimes
std::vector<Eigen::Matrix3d> Calculate_Lprimes(Lattice my_lattice) // Again may need to go with 
{
    Eigen::Matrix3d lattice=my_lattice.col_vector_matrix();
    int radius=1;
    std::vector<Eigen::Matrix3d> Lprimes;
    auto PS = calculate_gridpoints(lattice, radius);
    Eigen::Matrix3d MakeMatrix;
    for (auto p1 : PS)
    {
        for (auto p2 : PS)
        {
            for (auto p3 : PS)
            {
                MakeMatrix << p1(0), p2(0), p3(0), p1(1), p2(1), p3(1), p1(2), p2(2),
                    p3(2); // I may have these flipped. If Im supposed to have [p1][p2] then this woulr be correct
                Lprimes.push_back(MakeMatrix.inverse());
            }
        }
    }
    return Lprimes;
}

// Determine whether or not the calculated Symmetry operation is valid through comparison of S^T*S=I
bool is_symop_valid(Eigen::Matrix3d SymMatrix)
{
    auto Matrixcheck = SymMatrix.transpose() * SymMatrix;
    if (!Matrixcheck.isIdentity(.0005))
        return false;
    else
        return true;
}
// This function calculates the symmetry operations that are valid for a given lattice
std::vector<SymOp> Calculate_point_group(Lattice my_lattice) // Is the type symops?
{
    Eigen::Matrix3d lattice=my_lattice.col_vector_matrix();
    int radius=1;
    Eigen::Matrix3d dummy_matrix;
    Eigen::Vector3d dummy_tau;
    SymOp test_symop(dummy_matrix, dummy_tau);
    std::vector<SymOp> validsymops;
    auto Lprimes = Calculate_Lprimes(lattice);
    Eigen::Matrix3d SymmetryOp;
    for (auto Lp : Lprimes)
    {
        SymmetryOp = lattice * Lp;
        if (is_symop_valid(SymmetryOp))
        {
            test_symop.my_matrix=SymmetryOp;
	    test_symop.my_trans<<0,0,0; //would I need to push this back as well?
	    validsymops.push_back(test_symop);
        }
    }
    return validsymops;
}

//    bool is_SymOp_Unitary(SymOp SymMatrix); //Look up C++ general determinant calcs

//
bool group_is_closed(std::vector<Eigen::Matrix3d> SymMatrix) //
{
    Eigen::Matrix3d GroupMultiplication;
    for (auto S1 : SymMatrix)
    {
        for (auto S2 : SymMatrix)
        {
            GroupMultiplication = S1 * S2;
            compare_mat compare_matrix(GroupMultiplication);

            if (std::find_if(SymMatrix.begin(), SymMatrix.end(), compare_matrix) != SymMatrix.end())
                return true;
            else
                return false;
        }
    }
}


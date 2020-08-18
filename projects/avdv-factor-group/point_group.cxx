#include "../eigen-git-mirror/Eigen/Dense"
#include "./math.hpp"
#include "./symop.hpp"
#include "./symgroup.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace std;


std::vector<std::vector<double>> create_grid_pts(const Eigen::Matrix3d L)
{
    // generate a grid of coordinates of points with radius n.
    // each coordinate is vector triple of doubles in a 2n+1 by 2n+1 eigen matrix
    double n = 2.0;
    int l = 2 * n + 1;
    std::vector<std::vector<double>> grid;
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < l; j++)
        {
            for (int k = 0; k < l; k++)
            {
                std::vector<double> temp;
                double ii = i;
                double jj = j;
                double kk = k;
                temp.push_back((ii - n) * L(0, 0) + (jj - n) * L(0, 1) + (kk - n) * L(0, 2));
                temp.push_back((ii - n) * L(1, 0) + (jj - n) * L(1, 1) + (kk - n) * L(1, 2));
                temp.push_back((ii - n) * L(2, 0) + (jj - n) * L(2, 1) + (kk - n) * L(2, 2));
                grid.push_back(temp);
            }
        }
    }

    return grid;
}

std::vector<Eigen::Matrix3d> calc_L_primes(const std::vector<std::vector<double>> grid)
{
    // calculate the L-primes, which are sets of three coordinates or lattice vectors representing
    // possible sets of transformed Lattice vectors

    std::vector<Eigen::Matrix3d> L_prime;
    for (const auto p1 : grid)
    {
        for (const auto p2 : grid)
        {
            for (const auto p3 : grid)
            {
                Eigen::Matrix3d temp;
                temp(0, 0) = p1[0];
                temp(1, 0) = p1[1];
                temp(2, 0) = p1[2];
                temp(0, 1) = p2[0];
                temp(1, 1) = p2[1];
                temp(2, 1) = p2[2];
                temp(0, 2) = p3[0];
                temp(1, 2) = p3[1];
                temp(2, 2) = p3[2];
                L_prime.push_back(temp);
            }
        }
    }
    return L_prime;
}

SymGroup<SymOp, CartesianBinaryComparator_f> generate_point_group(const Eigen::Matrix3d& L, double tol)
{
    std::vector<SymOp> pt_group_list;
    auto pts = create_grid_pts(L);
    auto L_primes = calc_L_primes(pts);
    for (Eigen::Matrix3d LP : L_primes)
    {
        Eigen::Matrix3d cart_mat=LP * L.inverse();
        // TODO: Pass tolerance, don't hard code
        if (!is_unitary(cart_mat, tol))
        {
            continue;
        }

        pt_group_list.emplace_back(cart_mat);
    }
    CartesianBinaryComparator_f compare(tol);
    SymGroup<SymOp, CartesianBinaryComparator_f> pt_group(pt_group_list,compare);
    return pt_group;
}

SymGroup<SymOp, CartesianBinaryComparator_f> generate_point_group(const Lattice& lat, double tol)
{
    const Eigen::Matrix3d& column_matrix=lat.col_vector_matrix();
    return generate_point_group(column_matrix,tol);
}


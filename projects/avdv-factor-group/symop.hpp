#ifndef SYMOP_H
#define SYMOP_H

#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include "../../submodules/eigen-git-mirror/Eigen/Dense"
#include <string>
#include <vector>
#include "lattice.hpp"

class SymOp
{
public:
    SymOp(const Eigen::Matrix3d& cart_matrix, const Eigen::Vector3d& translation={0, 0, 0});
    Eigen::Vector3d get_translation() const;
    Eigen::Matrix3d get_cart_matrix() const;

private:
    Eigen::Matrix3d m_cart_matrix;
    Eigen::Vector3d m_translation;
};

SymOp operator*(const SymOp& lhs, const SymOp& rhs);


//TODO: Refactor. Just elimitate this, we only want the binary comparisons
//If you want unary comparisons, make a fancy template
//Change later
struct SymOpCompare_f
{

    SymOpCompare_f(SymOp input1, double tol=1e-5);
    bool operator()(const SymOp& element2) const;

private:
    const SymOp element1;
    double tol;
};

//TODO: Better name. This doesn't describe anything except the function signature
class CartesianBinaryComparator_f
{
    public:
            CartesianBinaryComparator_f(double tol);
            bool operator()(const SymOp& element1, const SymOp& element2) const;
    private:
            double tol;
};

class BinarySymOpPeriodicCompare_f
{
	public:
		BinarySymOpPeriodicCompare_f(const Lattice& lattice, double tol);
		bool operator()(const SymOp& element1, const SymOp& element2) const;
	private:
		double tol;
		Lattice m_lattice;
};

class BinarySymOpPeriodicMultiplier_f
{
    public:
        BinarySymOpPeriodicMultiplier_f(const Lattice& lattice, double tol);
        SymOp operator()(const SymOp& operation1, const SymOp& operation2) const;
    private:
        Lattice m_lattice;
        double tol;
};

bool operator==(const SymOp& lhs, const SymOp& rhs);

#endif

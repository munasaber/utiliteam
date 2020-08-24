#ifndef COORDINATE_H
#define COORDINATE_H
#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include "../../submodules/eigen-git-mirror/Eigen/Dense"
#include "lattice.hpp"
#include <string>
#include <vector>
#include "symop.hpp"
struct VectorCompare_f
{
	VectorCompare_f(const Eigen::Vector3d& vector, double prec);
	bool operator()(const Eigen::Vector3d& other) const;
	private:
	Eigen::Vector3d m_vector;
	double m_precision;
};

struct VectorPeriodicCompare_f
{
    VectorPeriodicCompare_f(const Eigen::Vector3d& vector, double prec, const Lattice& unit_cell);
    bool operator()(const Eigen::Vector3d& other) const;
private:
    Eigen::Vector3d m_vector;
    double m_precision;
    Lattice m_lattice;
};

Eigen::Vector3d operator*(const SymOp& transformation, const Eigen::Vector3d& vector);
#endif

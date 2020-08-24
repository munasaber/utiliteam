#include "vectorfunctions.hpp"
#include "lattice.hpp"

// Vector Compare functor
VectorCompare_f::VectorCompare_f(const Eigen::Vector3d& vector, double prec) : m_vector(vector), m_precision(prec) {}
bool VectorCompare_f::operator()(const Eigen::Vector3d& other) const
{
    if (m_vector.isApprox(other, m_precision))
    {
            return true;
    }
    return false;
}


VectorPeriodicCompare_f::VectorPeriodicCompare_f(const Eigen::Vector3d& vector, double prec, const Lattice& unit_cell) : m_vector(vector), m_precision(prec), m_lattice(unit_cell){
} 
bool VectorPeriodicCompare_f::operator()(const Eigen::Vector3d& other) const 
{	
    Eigen::Vector3d vector1 = convert_to_fractional(m_lattice, m_vector);
    Eigen::Vector3d vector2 = convert_to_fractional(m_lattice, other);
    Eigen::Vector3d distance_vector=vector1-vector2;
    
    for (int i=0; i<distance_vector.size();i++){
        distance_vector(i)=distance_vector(i)-std::round(distance_vector(i));
    }
    Eigen::Vector3d cartesian_distance_vector = convert_to_cartesian(m_lattice, distance_vector);
    return std::abs(cartesian_distance_vector.norm())<m_precision;
}

Eigen::Vector3d operator*(const SymOp& transformation, const Eigen::Vector3d& vector)
{
    Eigen::Matrix3d transformed_matrix= (transformation.get_cart_matrix());
    Eigen::Vector3d transformed = (transformed_matrix)*vector+transformation.get_translation();
    return transformed;
}


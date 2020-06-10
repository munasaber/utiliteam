#include "coordinate.hpp"
#include "lattice.hpp"
// Coordinate defintions
Coordinate::Coordinate(const Eigen::Vector3d& coord) : m_coord(coord) {}
Eigen::Vector3d Coordinate::get_cart_coordinate() const { return this->m_coord; }
double Coordinate::get_x() const { return this->m_coord(0); }
double Coordinate::get_y() const { return this->m_coord(1); }
double Coordinate::get_z() const { return this->m_coord(2); }
Eigen::Vector3d Coordinate::get_frac_coordinate(const Lattice& lattice) const{   
    return convert_to_fractional(lattice, m_coord);
} 

void Coordinate::coord_bring_within(const Lattice& lattice, double prec)
{
    
//    Eigen::Vector3d frac_coords;
    Eigen::Vector3d cart_coords = this->get_cart_coordinate();
    this->m_coord = bring_within(lattice, prec, cart_coords);
    //frac_coords = lattice.col_vector_matrix().inverse() * get_cart_coordinate(lattice);
/*    for (int i = 0; i < 3; ++i)
    {
        if (frac_coords(i) < -prec || frac_coords(i) >= 1-prec)		
        {
            frac_coords(i) = frac_coords(i) - floor(frac_coords(i)+prec);
        }
    }
    this->m_coord = lattice.col_vector_matrix() * frac_coords;*/
}
//Is this an altright location for this? I can't really think of another.


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
    Eigen::Vector3d transformed = (transformation.get_cart_matrix()) * (vector)+transformation.get_translation();
    return transformed;
}

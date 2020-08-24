#include "site.hpp"
#include <cmath>

// Defines position and type of atom in a crystal
Site::Site(const std::string atom_name, const Eigen::Vector3d& init_coord) : m_coord(init_coord), m_atom(atom_name) {}
std::string Site::get_atom() const { return this->m_atom; }
Eigen::Vector3d Site::get_eigen_coordinate() const { return this->m_coord; }
//Coordinate Site::get_coordinate() const { return this->m_coord; }

// Site Compare functor
SiteCompare_f::SiteCompare_f(const Site& site, double prec) : m_site(site), m_precision(prec) {}
bool SiteCompare_f::operator()(const Site& other) const
{
    if (m_site.get_eigen_coordinate().isApprox(other.get_eigen_coordinate(), m_precision))
    {
        if (m_site.get_atom() == other.get_atom())
        {
            return true;
        }
    }
    return false;
}



//Site compare taking into account perioduc boundary conditions
SitePeriodicCompare_f::SitePeriodicCompare_f(Site site, double prec, const Lattice& unit_cell) : m_site(site), m_precision(prec), m_lattice(unit_cell){
} 
bool SitePeriodicCompare_f::operator()(Site other) 
{	
    VectorPeriodicCompare_f compare(m_site.get_eigen_coordinate(), m_precision, m_lattice);	
    return m_site.get_atom()==other.get_atom() && compare(other.get_eigen_coordinate());
}


Site operator*(const SymOp& transformation, const Site& site)
{
    Eigen::Vector3d site_coords=site.get_eigen_coordinate();
    Eigen::Vector3d transformed = transformation*(site_coords)+transformation.get_translation();
    return Site(site.get_atom(), transformed);
}

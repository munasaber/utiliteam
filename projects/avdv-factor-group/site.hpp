#ifndef SITE_H
#define SITE_H
#include "vectorfunctions.hpp"
#include "lattice.hpp"
#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include "../../submodules/eigen-git-mirror/Eigen/Dense"
#include <string>
#include <vector>
#include "symop.hpp"
class Site
{
public:
    Site(const std::string atom_name, const Eigen::Vector3d& init_coord);
    std::string get_atom() const;
    Eigen::Vector3d get_eigen_coordinate() const;
    Eigen::Vector3d m_coord;
private:
    std::string m_atom;
   // Coordinate m_coord;
};


//TODO: make binary
/// Compares two sites to see if coordinates and atom types match
struct SiteCompare_f
{
    SiteCompare_f(const Site& site, double prec);
    bool operator()(const Site& other) const;

private:
    Site m_site;
    double m_precision;
};

/// Compares two sites and returns true if the types match and
/// the coordinates are the same after applying unit cell tranlsations

struct SitePeriodicCompare_f
{
    SitePeriodicCompare_f(Site site, double prec, const Lattice& unit_cell);
    bool operator()(Site other);
private:
    Site m_site;
    double m_precision;
    Lattice m_lattice;
};

Site operator*(const SymOp& transformation, const Site& site);
#endif

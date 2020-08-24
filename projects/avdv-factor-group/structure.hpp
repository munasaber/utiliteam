#ifndef STRUCTURE_H
#define STRUCTURE_H 
#include "site.hpp"
#include "lattice.hpp"
#include "vectorfunctions.hpp"
#include "symop.hpp"
#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include "../../submodules/eigen-git-mirror/Eigen/Dense"
#include <string>
#include <vector>
class Structure
{
public:
    Structure(const Lattice& lattice, const std::vector<Site>& sites);
    Lattice get_lattice() const;
    std::vector<Site> get_sites() const;

private:
    Lattice m_lattice;
    std::vector<Site> m_sites;
};
#endif

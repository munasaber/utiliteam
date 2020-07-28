#include "write_to_poscar.hpp"
#include <unordered_map>
void write_to_poscar(const Structure& structure, std::ofstream& poscar_stream)
{	
	poscar_stream<<"POSCAR"<<'\n';
	poscar_stream<<1.0<<'\n';
	poscar_stream<<structure.get_lattice().lattice_vector(0).transpose()<<'\n'<<structure.get_lattice().lattice_vector(1).transpose()<<'\n'<<structure.get_lattice().lattice_vector(2).transpose()<<'\n';
 	std::unordered_map<std::string, std::vector<Eigen::Vector3d>> atom_name_map;
	for (const auto& site: structure.get_sites())
	{
		Eigen::Vector3d vect=site.get_eigen_coordinate();
		std::string name=site.get_atom();
		atom_name_map[name].push_back(vect);
	}
	for (const auto& atom:atom_name_map)
	{
		poscar_stream<<atom.first<<' ';
	}
	poscar_stream<<std::endl;	
	for (const auto& atom:atom_name_map)
	{
		poscar_stream<<atom.second.size()<<' ';
	}	
	poscar_stream<<std::endl;
	poscar_stream<<"Direct"<<'\n';//change to sifted_grid
	for (const auto& atom: atom_name_map)
	{
		for (const auto& atom_vector: atom.second)
		{
			poscar_stream<<convert_to_fractional(structure.get_lattice(), atom_vector).transpose()<<std::endl;
		}
	}
	poscar_stream.close();
}

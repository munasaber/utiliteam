#include "structure.hpp"
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

Structure read_poscar(const std::string &poscar_path) {
  std::ifstream file;
  file.open(poscar_path, std::ios::in);
  std::string line;

  // Comment line
  std::getline(file, line);

  // Scaling line
  std::getline(file, line);
  std::istringstream ss0(line);
  double scaling;
  ss0 >> scaling;

  // Eigen Matrix line
  Eigen::Matrix3d lat_row_matrix;
  for (int i = 0; i < 3; ++i) {
    std::getline(file, line);
    std::istringstream ss1(line);
    ss1 >> lat_row_matrix(i, 0) >> lat_row_matrix(i, 1) >> lat_row_matrix(i, 2);
  }

  // species line
  std::getline(file, line);
  std::vector<std::string> species;
  std::istringstream ss2(line);
  std::string specie;
  while (ss2 >> specie) {
    species.push_back(specie);
  }

  // Number of species line
  std::getline(file, line);
  std::vector<int> NumSpecies;
  std::istringstream ss3(line);
  int spec;
  while (ss3 >> spec) {
    NumSpecies.push_back(spec);
  }

  // Coord-type line
  std::getline(file, line);
  std::string coord_type;
  std::istringstream ss4(line);
  ss4 >> coord_type;
  

  // Raw coordinates
  Eigen::Vector3d coord;
  std::vector<Eigen::Vector3d> raw_coordinate_values;
  char m= std::tolower(coord_type[0]);
  if (m=='c')
  {	  
  	while (std::getline(file, line)) 
	{
    		std::istringstream ss5(line);
    		ss5 >> coord(0) >> coord(1) >> coord(2);
    		raw_coordinate_values.push_back(coord);
  	}
  }
  //else if (m=='d')
  if (m=='d')
  {
  	while (std::getline(file, line)) 
	{
    		std::istringstream ss5(line);
    		ss5 >> coord(0) >> coord(1) >> coord(2);
		Eigen::Vector3d cart_coord=lat_row_matrix.transpose()*coord;
    		raw_coordinate_values.push_back(cart_coord);
  	}
  }

  //if coordtype='Cartesian' convert to cart
  

  // make Lattice
  Lattice latt(lat_row_matrix.row(0),lat_row_matrix.row(1),lat_row_matrix.row(2));
  // making the sites, creating one site per atom
  std::vector<Site> Sites;
  int t = 0;
  for (int i = 0; i < species.size(); i++)
    for (int j = 0; j < NumSpecies[i]; j++) {
      Eigen::Vector3d temp(raw_coordinate_values[t]);
      Site Single_site(species[i], temp);
      Sites.push_back(Single_site);
      t = t + 1;
    }

  // making the structure
  // for instance printing first row
  return Structure(latt, Sites);
};


#g++ interstitials.cpp xtal_classes.cpp symmetry_operations.cpp tests.cpp -o run
g++ --std=c++17 ../avdv-factor-group/lattice.cpp ../avdv-factor-group/site.cpp ../avdv-factor-group/coordinate.cpp ../avdv-factor-group/symop.cpp ../avdv-factor-group/structure.cpp ../avdv-factor-group/tests.cpp ../avdv-factor-group/symgroup.hpp ../avdv-factor-group/factor_group.cpp ../avdv-factor-group/point_group.cxx interstitial_mesh.cpp interstitial_mesh_test.cpp -o interstitial_test
g++ --std=c++17 -g -Og ../avdv-factor-group/lattice.cpp ../avdv-factor-group/site.cpp ../avdv-factor-group/coordinate.cpp ../avdv-factor-group/symop.cpp ../avdv-factor-group/structure.cpp write_to_poscar.cpp ../avdv-factor-group/symgroup.hpp ../avdv-factor-group/factor_group.cpp ../avdv-factor-group/point_group.cxx interstitial_mesh.cpp main_pnb9o25.cpp -o test_main_pnb9o25
#g++ --std=c++17 -g -Og ../avdv-factor-group/lattice.cpp ../avdv-factor-group/site.cpp ../avdv-factor-group/coordinate.cpp ../avdv-factor-group/symop.cpp ../avdv-factor-group/structure.cpp write_to_poscar.cpp ../avdv-factor-group/symgroup.hpp ../avdv-factor-group/factor_group.cpp ../avdv-factor-group/point_group.cxx interstitial_mesh.cpp main_fcc.cpp -o test_main_fcc
#g++ --std=c++17 -g -Og ../avdv-factor-group/lattice.cpp ../avdv-factor-group/site.cpp ../avdv-factor-group/coordinate.cpp ../avdv-factor-group/symop.cpp ../avdv-factor-group/structure.cpp write_to_poscar.cpp ../avdv-factor-group/symgroup.hpp ../avdv-factor-group/factor_group.cpp ../avdv-factor-group/point_group.cxx interstitial_mesh.cpp main_cubic.cpp -o test_main_cubic

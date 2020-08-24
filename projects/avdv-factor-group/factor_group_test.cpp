#include "./factor_group.hpp"
#include "./symgroup.hpp"
#include "./symop.hpp"
#include "./io.hpp"
#include "./tests.hpp"

//Lattice identity_lattice()
//{
//    //TODO
//}

bool test_transform_basis_inversion()
{
    SymOp inversion(Eigen::Matrix3d::Identity() * -1);
    std::vector<Site> basis{Site("Na", Eigen::Vector3d(1, 1, 1))};

    auto inverted_basis = transform_basis(inversion, basis);
    SiteCompare_f compare(inverted_basis.at(0), 1e-5);

    return compare(Site("Na", Eigen::Vector3d(-1, -1, -1)));
}

bool test_transform_basis_glide()
{
    Eigen::Matrix3d mirror_mat = Eigen::Matrix3d::Identity();
    mirror_mat(2, 2) = -1;
    SymOp glide(mirror_mat, Eigen::Vector3d(0.2, 0.2, 0));

    std::vector<Site> basis{Site("Na", Eigen::Vector3d(1, 1, 1))};
    auto glided_basis = transform_basis(glide, basis);
    SiteCompare_f compare(glided_basis.at(0), 1e-5);

    return compare(Site("Na", Eigen::Vector3d(1.2, 1.2, -1)));
}

bool test_transform_basis_size()
{
    Site s("Na", Eigen::Vector3d(1, 1, 1));
    std::vector<Site> basis(100, s);
    SymOp identity(Eigen::Matrix3d::Identity());
    auto same_basis = transform_basis(identity, basis);
    return same_basis.size() == basis.size();
}

bool basis_maps_onto_itself_test(double tol)
{
    Eigen::Matrix3d identity=Eigen::Matrix3d::Identity();
    Lattice lat(identity.col(0),identity.col(1),identity.col(2));
    Site na1("Na",Eigen::Vector3d::Zero());
    Site na2("Na",0.5*Eigen::Vector3d::Ones());
    Site li1("Li",-0.25*Eigen::Vector3d::Ones());
    std::vector<Site> basis{na1,na2,li1};
    std::vector<Site> basis_shuffle{li1,na2,na1};

    return basis_maps_onto_itself(basis,basis,lat, tol) &&
    basis_maps_onto_itself(basis,basis_shuffle,lat, tol);
    
}

bool basis_maps_onto_itself_periodically(double tol)
{
    Eigen::Matrix3d identity=Eigen::Matrix3d::Identity();
    Lattice lat(identity.col(0),identity.col(1),identity.col(2));
    Site na1("Na",Eigen::Vector3d::Zero());
    Site na2("Na",0.5*Eigen::Vector3d::Ones());
    Site li1("Li",-0.25*Eigen::Vector3d::Ones());

    std::vector<Site> basis{na1,na2,li1};
    Site periodic_na2("Na",1.5*Eigen::Vector3d::Ones());
    std::vector<Site> basis_periodic{na1,periodic_na2,li1};

    return basis_maps_onto_itself(basis,basis_periodic,lat, tol);
}

bool basis_doesnt_map_onto_itself(double tol)
{
    Eigen::Matrix3d identity=Eigen::Matrix3d::Identity();
    Lattice lat(identity.col(0),identity.col(1),identity.col(2));
    Site na1("Na",Eigen::Vector3d::Zero());
    Site na2("Na",0.5*Eigen::Vector3d::Ones());
    Site li1("Li",-0.25*Eigen::Vector3d::Ones());
    Site na3("Na",-0.25*Eigen::Vector3d::Ones());

    return !basis_maps_onto_itself({na1,na2,li1},{na1,na2,na3},lat, tol);
}

bool basis_doesnt_map_onto_itself_size_mismatch(double tol)
{
    Eigen::Matrix3d identity=Eigen::Matrix3d::Identity();
    Lattice lat(identity.col(0),identity.col(1),identity.col(2));
    Site na1("Na",(Eigen::Vector3d::Zero()));
    Site na2("Na",(0.5*Eigen::Vector3d::Ones()));
    Site li1("Li",(-0.25*Eigen::Vector3d::Ones()));
    Site na3("Na",(-0.25*Eigen::Vector3d::Ones()));

    return !basis_maps_onto_itself({na1,na2,li1},{na1},lat, PREC) && !basis_maps_onto_itself({na1},{na1,na2,li1},lat, tol);
}

bool generate_translations_works(double tol)
{
	Site original_basis("Na", (Eigen::Vector3d(0.0, 0.0, 0.0)));
	Eigen::Vector3d expected_translations1(Eigen::Vector3d(0.2, 0.2, 0.2));
	Eigen::Vector3d expected_translations2(Eigen::Vector3d(0.4, 0.2, 0.2));
	Eigen::Vector3d expected_translations3(Eigen::Vector3d(0.6, 0.2, 0.2));
        Site transformed_basis_site1("Na", expected_translations1);
        Site transformed_basis_site2("Na", expected_translations2);
        Site transformed_basis_site3("Na", expected_translations3);
	std::vector<Site> total_transformed{transformed_basis_site1, transformed_basis_site2, transformed_basis_site3};
	std::vector<Eigen::Vector3d> my_translations=generate_translations(original_basis, total_transformed);
	std::vector<Eigen::Vector3d> vector_expected_translations{expected_translations1*-1, expected_translations2*-1, expected_translations3*-1};
	if (my_translations.size()!=3)
	{
		std::cout<<"wrong size"<<std::endl;
		std::cout<<my_translations.size()<<std::endl;
		return false;
	}
/*	std::cout<<"my_translations";
	for(const auto& translation: my_translations)
	{
	
		std::cout<<translation<<std::endl;
	}
	std::cout<<"expected_translation";*/

	for(const auto& translation: vector_expected_translations)
	{
	
//		std::cout<<translation<<std::endl;
	}

	for(const auto& translation: my_translations)
	{
		auto compare_translations=[translation, tol](Eigen::Vector3d expected_translation)	
		{
			return translation.isApprox(expected_translation, tol);
		};
		 if (find_if(vector_expected_translations.begin(), vector_expected_translations.end(), compare_translations)==vector_expected_translations.end())
			{	
		 	std::cout<<"generated translations not within expected_translations";
			std::cout<<translation<<std::endl;
			return false;
			}
	}

	return true;
}

bool test_fcc_factor_group(double tol)
{
    Structure fcc=read_poscar("./test_files/fcc.vasp");
    auto fg=generate_factor_group(fcc, tol);
    
    return fg.operations().size()==48;
}

bool test_hcp_factor_group(double tol)
{
    Structure hcp=read_poscar("./test_files/mg.vasp");
    auto fg=generate_factor_group(hcp, tol);
    
    return fg.operations().size()==24;
}

bool test_FeLi2Se2_factor_group(double tol)
{
    Structure fcc=read_poscar("./test_files/FeLi2Se2.vasp");
    auto fg=generate_factor_group(fcc, tol); 
    std::cout<<"DEBUGGING: feli2se2 fg.size() is "<<fg.operations().size()<<std::endl;
  //std::cout<<"Hellooooooooooo\n";
    for(const SymOp& op : fg.operations())
    {
        std::cout<<op.get_cart_matrix()<<"\n\n";
        std::cout<<op.get_translation().transpose()<<"\n\n";
    }
    return fg.operations().size()==6;
}

bool test_zro2_factor_group(double tol)
{
    Structure zro2=read_poscar("./test_files/POSCAR.ZrO2");
    auto fg=generate_factor_group(zro2, tol); 
    std::cout<<"DEBUGGING: ZrO2 fg.size() is "<<fg.operations().size()<<std::endl;
    for(const SymOp& op : fg.operations())
    {
        std::cout<<op.get_cart_matrix()<<"\n\n";
        std::cout<<op.get_translation().transpose()<<"\n\n";
    }
    return fg.operations().size()==4;
}
bool test_pnb9o25_factor_group(double tol)
{
    Structure pnb9o25=read_poscar("./test_files/pnb9o25.vasp");
    auto fg=generate_factor_group(pnb9o25, tol);
    std::cout<<"DEBUGGING: pnb9o25 fg.size() is "<<fg.operations().size()<<std::endl;
//    std::cout<<"Hellooooooooooo\n";

    //TODO: Make SymGroup work for range based loops (you'll have to forward begin/end iterators)
    for(const SymOp& op : fg.operations())
    {
        std::cout<<op.get_cart_matrix()<<"\n\n";
        std::cout<<op.get_translation().transpose()<<"\n\n";
    }

    return fg.operations().size()==4;
}

int main()
{
    std::cout<<"---- Running Factor Group Tests ----"<<std::endl;
    std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
    double tol= 1e-4;
    EXPECT_TRUE(test_transform_basis_inversion(), "Transform basis inversion");
    EXPECT_TRUE(test_transform_basis_glide(), "Transform basis glide");
    EXPECT_TRUE(test_transform_basis_size(), "Transform basis size");

    EXPECT_TRUE(basis_maps_onto_itself_test(tol), "Basis maps");
    EXPECT_TRUE(basis_maps_onto_itself_periodically(tol), "Basis maps periodic");
    EXPECT_TRUE(basis_doesnt_map_onto_itself(tol), "Basis shouldnt map, has different atom types");
    EXPECT_TRUE(generate_translations_works(tol), "Generate translations should generate translations");
    EXPECT_TRUE(test_fcc_factor_group(tol), "FCC factor group has 48 operations");
    EXPECT_TRUE(test_hcp_factor_group(tol), "hcp factor group has 24 operations");
    EXPECT_TRUE(test_pnb9o25_factor_group(tol), "pnb9o25 factor group has 4 operations");
//    EXPECT_TRUE(test_FeLi2Se2_factor_group(tol), "FeLi2Se2 factor group has 48 operations");
    EXPECT_TRUE(test_zro2_factor_group(tol), "ZrO2 factor group has 4 operations");
    return 0;
}

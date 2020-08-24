#include <algorithm>
#include <vector>
//#include "./symop.hpp"
//#include "./site.hpp"
//#include "./symgroup.hpp"
//#include "./structure.hpp"
#include "./factor_group.hpp" 

std::vector<Site> transform_basis(const SymOp& symop, const std::vector<Site>& basis)
{
    std::vector<Site> transformed_basis;
    for(const Site& s : basis)
    {
        transformed_basis.emplace_back(symop*s);
    }
    return transformed_basis;
}

//TODO: tolerance parameter makes it not match delclaration anymore
bool basis_maps_onto_itself(const std::vector<Site>& original_basis, const std::vector<Site>& transformed_basis, const Lattice& lattice, double tol)
{
    if(original_basis.size()!=transformed_basis.size())
    {
	 return false;
	
    }

    for(const Site& transformed_s : transformed_basis)
    {
        SitePeriodicCompare_f compare_to_transformed(transformed_s,tol,lattice);
        
	if(std::find_if(original_basis.begin(),original_basis.end(),compare_to_transformed)==original_basis.end())
        {
	return false;
	    
        }

    }

    return true;
}


std::vector<Eigen::Vector3d> generate_translations(const Site& original_basis_site, const std::vector<Site>& transformed_basis)
{
	std::vector<Eigen::Vector3d> total_trans;
	for (const Site& transformed_basis_site: transformed_basis)
	{
		Eigen::Vector3d trans= original_basis_site.get_eigen_coordinate()-transformed_basis_site.get_eigen_coordinate();
		total_trans.push_back(trans);
	}
	return total_trans;
}



SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f> generate_factor_group(const Structure& struc, double tol)
{
    SymGroup<SymOp, CartesianBinaryComparator_f>  point_group=generate_point_group(struc.get_lattice(), tol);
    const auto& basis=struc.get_sites(); 

    //make empty sym group
    BinarySymOpPeriodicCompare_f comparison(struc.get_lattice(), tol);
    BinarySymOpPeriodicMultiplier_f mult_op(struc.get_lattice(), tol);
    SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f> factor_group({},comparison, mult_op);
    for(const SymOp& point_op : point_group.operations())
   {	
        auto transformed_basis=transform_basis(point_op,basis);
        std::vector<Eigen::Vector3d> all_translations=generate_translations(basis[0],transformed_basis);

        for(const Eigen::Vector3d& translation : all_translations)
        {
            SymOp symop_translation(Eigen::Matrix3d::Identity(),translation);
            auto transformed_translated_basis=transform_basis(symop_translation,transformed_basis);
            if(basis_maps_onto_itself(basis,transformed_translated_basis,struc.get_lattice(),tol))
            {
	        	SymOp factor_group_symop(point_op.get_cart_matrix(), translation);
                factor_group.insert(factor_group_symop);
            }
        }
    }

    return factor_group;
}

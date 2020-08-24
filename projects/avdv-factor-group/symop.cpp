#include "symop.hpp"
#include "site.hpp"
#include  "vectorfunctions.hpp"
#include "lattice.hpp"

SymOp::SymOp(const Eigen::Matrix3d& cart_matrix, const Eigen::Vector3d& translation)
    : m_cart_matrix(cart_matrix), m_translation(translation)
{
}
Eigen::Vector3d SymOp::get_translation() const { return this->m_translation; }
Eigen::Matrix3d SymOp::get_cart_matrix() const { return this->m_cart_matrix; }

SymOp operator*(const SymOp& lhs, const SymOp& rhs)
{
    Eigen::Vector3d translation = lhs.get_cart_matrix() * rhs.get_translation() + lhs.get_translation();
    Eigen::Matrix3d product = lhs.get_cart_matrix() * rhs.get_cart_matrix();
    SymOp symop_product(product, translation);
    return symop_product;
}


//TODO: Make this disappear
SymOpCompare_f::SymOpCompare_f(SymOp input1, double tol) : element1(input1), tol(tol) {}

bool SymOpCompare_f::operator()(const SymOp& element2) const
{
    return (element1.get_cart_matrix().isApprox(element2.get_cart_matrix(), tol) &&
            (element1.get_translation().isApprox(element2.get_translation(), tol)));
}

CartesianBinaryComparator_f::CartesianBinaryComparator_f(double tol) : tol(tol) {}

bool CartesianBinaryComparator_f::operator()(const SymOp& element1, const SymOp& element2) const
{

   return (element1.get_cart_matrix().isApprox(element2.get_cart_matrix(), tol) &&
            (element1.get_translation().isApprox(element2.get_translation(), tol)));
}


BinarySymOpPeriodicCompare_f::BinarySymOpPeriodicCompare_f(const Lattice& lattice, double tol) : m_lattice(lattice), tol(tol) {}

bool BinarySymOpPeriodicCompare_f::operator()(const SymOp& element1, const SymOp& element2) const
{    
	Site temp_site1 = Site(std::string("xx"), element1.get_translation());
	Site temp_site2 = Site(std::string("xx"), element2.get_translation());
    SitePeriodicCompare_f translation_comparison(temp_site1, tol, m_lattice);
    
    SymOp symop1(element1.get_cart_matrix());
	SymOp symop2(element2.get_cart_matrix());
	CartesianBinaryComparator_f compare(tol);

	return compare(symop1, symop2) && translation_comparison(temp_site2);
}

BinarySymOpPeriodicMultiplier_f::BinarySymOpPeriodicMultiplier_f(const Lattice& lattice, double tol) : m_lattice(lattice), tol(tol) {}

SymOp BinarySymOpPeriodicMultiplier_f::operator()(const SymOp& operation1, const SymOp& operation2) const 
{
    SymOp full_operation_product = operation1 * operation2;
    Eigen::Vector3d op_product_periodic_tranlation = bring_within(m_lattice, tol, full_operation_product.get_translation());
    SymOp final_product(full_operation_product.get_cart_matrix(), op_product_periodic_tranlation);
    return final_product;
}

bool operator==(const SymOp& lhs, const SymOp& rhs)
{
    CartesianBinaryComparator_f binarycompare(1e-6);
    return binarycompare(lhs, rhs);
}


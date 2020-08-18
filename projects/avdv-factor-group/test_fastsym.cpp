#include "./fastsymmetry.hpp"
#include "./symop.hpp"
#include "./point_group.hpp"
#include "./io.hpp"
#include "./structure.hpp"
#include <algorithm>
#include <iterator>
#include <string>
#include "./tests.hpp"
#include "./lattice.hpp"
#include "./site.hpp"
#include "./coordinate.hpp"

#define PREC 1e-6

void EXPECT_EQUAL_INT(int lhs, int  rhs, std::string test_name)
{ if(lhs!=rhs){ std::cout<<"FAILED TEST: "<<test_name<<std::endl;}
  if(lhs==rhs){std::cout<<"PASSED TEST: "<<test_name<<std::endl;}
}

int main(int argc, char *argv[])
{
    std::cout<<"---- Runnning FastSymmetry Tests ----"<<std::endl;
    std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
    
    if(argc < 2)
    {
        std::cout<<"Need a poscar file"<<std::endl;
        return 1;
    }
    
    //TODO: Use consistent style for naming classes (Youve been using CamelCase so far)
    Structure structure=read_poscar(argv[1]);
    SymGroup<SymOp, CartesianBinaryComparator_f> pt_group = generate_point_group(structure.get_lattice(), PREC);
    const auto& pt_group_operations=pt_group.operations();
    MultTable multiplication_table = make_multiplication_table(pt_group_operations, PREC);
    int group_sz=pt_group.operations().size();
//etst multiplication table construction
    EXPECT_EQUAL_INT(group_sz, multiplication_table.size(), "Multiplcation table row size check");
    EXPECT_EQUAL_INT(group_sz, multiplication_table[0].size(), "Multiplcation table column size check");
    EXPECT_EQUAL_INT( multiplication_table[0].size(),  multiplication_table.size(),"Multiplcation table is square check");

//test abstract Operation construction
    std::shared_ptr<MultTable> pt_group_multiplication_table_ptr = std::make_shared<MultTable>(multiplication_table);
    AbstractSymOp abstract_op1(1, pt_group_multiplication_table_ptr);
    EXPECT_TRUE(1==abstract_op1.get_id() && pt_group_multiplication_table_ptr==abstract_op1.mult_table_ptr(), "Abstract Operation Construction Test");

//test multiplication
    AbstractSymOp abstract_op4(4, pt_group_multiplication_table_ptr);
    AbstractSymOp product_op=abstract_op1*abstract_op4;
    EXPECT_EQUAL_INT(product_op.get_id(), multiplication_table[1][4],"Abstract Operation Multiplication Test");

//test compare abstract Ops
    AbstractSymOpCompare_f compare(abstract_op1);
    EXPECT_TRUE(!compare(abstract_op4), "Compared two different abstract ops");
    EXPECT_TRUE(compare(abstract_op1), "Compared two identical abstract ops");

//test transformation
    SymGroup<AbstractSymOp, BinaryAbstractComparator_f> abstract_pt_group = transform_representation(pt_group, PREC);
    EXPECT_TRUE(pt_group.operations().size()== abstract_pt_group.operations().size(), "Abstract Group and original group are same size");
    EXPECT_EQUAL_INT(1, abstract_pt_group.operations()[1].get_id(), "Check abstract id is index");
    EXPECT_TRUE(abstract_pt_group.operations()[1].mult_table_ptr()==abstract_pt_group.operations()[0].mult_table_ptr(), "check multiplication table ptrs are same in different group elements");

    //TODO:
    //Run through every multiplication, and make sure the abstract multiplications are
    //consistent with the Cartesian ones
    for(int i=0; i<pt_group_operations.size(); ++i)
    {
        for(int j=0; j<pt_group_operations.size(); ++j)
        {
            SymOp product=pt_group_operations[i]*pt_group_operations[j];
            SymOpCompare_f matches_product(product, PREC);
            auto product_it=std::find_if(pt_group_operations.begin(), pt_group_operations.end(), matches_product);
            int product_ix=std::distance(pt_group_operations.begin(),product_it);

            EXPECT_EQUAL_INT(product_ix,multiplication_table[i][j],"Check if multiplication table is consistent");
        }
    }

    return 0;
}


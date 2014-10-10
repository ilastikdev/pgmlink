#define BOOST_TEST_MODULE constraint_test

#include <vector>
#include <map>
#include <set>
#include <stdio.h>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/bind.hpp>
#include <boost/serialization/serialization.hpp>

#include <opengm/inference/icm.hxx>
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>

#include "pgmlink/constraint_pool.hxx"

using namespace pgmlink::pgm;

BOOST_AUTO_TEST_CASE(IncomingFunctionTest)
{
    // T1, T2, V nodes with label space 3
    std::vector<size_t> shape = {3, 3, 3};

    IncomingConstraintFunction<double, size_t, size_t> constraint_func(shape.begin(), shape.end());
    constraint_func.set_forbidden_energy(200.0);
    std::vector<size_t> labeling = {1,2,3};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    labeling = {1,2,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 200.0);

    labeling = {1,0,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    labeling = {1,1,2};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    labeling = {1,1,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 200.0);
}

BOOST_AUTO_TEST_CASE(OutgoingFunctionTest)
{
    // A, D, T1, T2
    std::vector<size_t> shape = {3, 3, 3, 3};

    OutgoingConstraintFunction<double, size_t, size_t> constraint_func(shape.begin(), shape.end());
    constraint_func.set_forbidden_energy(200.0);
    std::vector<size_t> labeling = {3,0,2,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    labeling = {1,0,2,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 200.0);

    labeling = {1,1,1,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    labeling = {1,1,2,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    labeling = {2,1,1,2}; // division not allowed when A > 1!
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 200.0);

    // test parameter effect
    constraint_func.set_with_divisions(false);

    labeling = {1,1,1,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 200.0);
}

BOOST_AUTO_TEST_CASE(DetectionFunctionTest)
{
    // V,A nodes with label space 3
    std::vector<size_t> shape = {3, 3};

    DetectionConstraintFunction<double, size_t, size_t> constraint_func(shape.begin(), shape.end());
    constraint_func.set_forbidden_energy(200.0);
    std::vector<size_t> labeling = {1,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    labeling = {1,2};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 200.0);

    labeling = {1,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    labeling = {0,2};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    labeling = {2,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 200.0);

    labeling = {0,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    // test the parameter effect
    constraint_func.set_with_appearance(false);
    labeling = {0,2};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 200.0);

    labeling = {1,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    constraint_func.set_with_appearance(true);
    constraint_func.set_with_disappearance(false);
    labeling = {0,2};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);

    labeling = {1,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 200.0);

    constraint_func.set_with_disappearance(true);
    constraint_func.set_with_misdetections(false);
    labeling = {0,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 200.0);

    labeling = {1,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin()), 0.0);
}

BOOST_AUTO_TEST_CASE(ConstraintPool_Size_Test)
{
    std::vector<size_t> dummy_vars = {1,2,3};

    ConstraintPool cp;
    cp.add_constraint(ConstraintPool::IncomingConstraint(dummy_vars, 4));
    cp.add_constraint(ConstraintPool::IncomingConstraint(dummy_vars, 4));
    cp.add_constraint(ConstraintPool::IncomingConstraint(dummy_vars, 4));

    BOOST_CHECK_EQUAL(cp.get_num_constraints(), 3);

    cp.add_constraint(ConstraintPool::OutgoingConstraint(1,2, dummy_vars));
    cp.add_constraint(ConstraintPool::OutgoingConstraint(1,2, dummy_vars));

    BOOST_CHECK_EQUAL(cp.get_num_constraints(), 5);

    cp.add_constraint(ConstraintPool::DetectionConstraint(1,2));
    cp.add_constraint(ConstraintPool::DetectionConstraint(1,3));
    cp.add_constraint(ConstraintPool::DetectionConstraint(8,2));

    BOOST_CHECK_EQUAL(cp.get_num_constraints(), 8);
}

BOOST_AUTO_TEST_CASE(ConstraintPool_Incoming_Factor_Test)
{
    OpengmModelDeprecated::ogmGraphicalModel model;
    model.addVariable(3);
    model.addVariable(3);
    model.addVariable(3);

    std::vector<size_t> labeling = {1,1,1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0.0);

    ConstraintPool cp(200.0);
    std::vector<size_t> indices = {0,1};
    cp.add_constraint(ConstraintPool::IncomingConstraint(indices, 2));

    opengm::ICM<OpengmModelDeprecated::ogmGraphicalModel, OpengmModelDeprecated::ogmAccumulator> inf(model);
    cp.add_constraints_to_problem(model, inf);

    // test wrong labelings
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    labeling = {1, 2, 2};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    // test good labelings
    labeling = {1, 1, 2};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);

    labeling = {0, 1, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);
}

BOOST_AUTO_TEST_CASE(ConstraintPool_Outgoing_Factor_Test)
{
    OpengmModelDeprecated::ogmGraphicalModel model;
    model.addVariable(3); // A
    model.addVariable(2); // D
    model.addVariable(3); // T
    model.addVariable(3); // T

    std::vector<size_t> labeling = {1,1,1,1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0.0);

    ConstraintPool cp(200.0);
    std::vector<size_t> indices = {2,3};
    cp.add_constraint(ConstraintPool::OutgoingConstraint(0, 1, indices));

    opengm::ICM<OpengmModelDeprecated::ogmGraphicalModel, OpengmModelDeprecated::ogmAccumulator> inf(model);
    cp.add_constraints_to_problem(model, inf);

    // test wrong labelings
    labeling = {2, 1, 2, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    labeling = {0, 1, 1, 0};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    labeling = {1, 1, 2, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    // test good labelings
    labeling = {1, 1, 2, 0};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);

    labeling = {2, 0, 1, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);
}

BOOST_AUTO_TEST_CASE(ConstraintPool_Outgoing_Factor_No_Division_Node_Test)
{
    OpengmModelDeprecated::ogmGraphicalModel model;
    model.addVariable(3); // A
    model.addVariable(3); // T
    model.addVariable(3); // T

    std::vector<size_t> labeling = {1,1,1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0.0);

    ConstraintPool cp(200.0);
    std::vector<size_t> indices = {1,2};
    cp.add_constraint(ConstraintPool::OutgoingConstraint(0, -1, indices));

    opengm::ICM<OpengmModelDeprecated::ogmGraphicalModel, OpengmModelDeprecated::ogmAccumulator> inf(model);
    cp.add_constraints_to_problem(model, inf);

    // test wrong labelings
    labeling = {2, 1, 2};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    labeling = {0, 1, 0};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    labeling = {1, 2, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    // test good labelings
    labeling = {1, 1, 0};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);

    labeling = {2, 1, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);
}

BOOST_AUTO_TEST_CASE(ConstraintPool_Detection_Factor_Test)
{
    OpengmModelDeprecated::ogmGraphicalModel model;
    model.addVariable(3); // A
    model.addVariable(3); // V

    std::vector<size_t> labeling = {1,2};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0.0);

    ConstraintPool cp(200.0);
    cp.add_constraint(ConstraintPool::DetectionConstraint(0, 1));

    opengm::ICM<OpengmModelDeprecated::ogmGraphicalModel, OpengmModelDeprecated::ogmAccumulator> inf(model);
    cp.add_constraints_to_problem(model, inf);

    // test wrong labelings
    labeling = {2, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    labeling = {1, 2};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    // test good labelings
    labeling = {1, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);

    labeling = {2, 0};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);

    labeling = {0, 2};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);
}

BOOST_AUTO_TEST_CASE(ConstraintPool_Serialization_Test)
{
    OpengmModelDeprecated::ogmGraphicalModel model;
    model.addVariable(3); // T 0
    model.addVariable(3); // T 1
    model.addVariable(3); // V 2
    model.addVariable(3); // A 3
    model.addVariable(2); // D 4
    model.addVariable(3); // T 5
    model.addVariable(3); // T 6

    // save model in this state:
    std::string filename(tmpnam(NULL));
    opengm::hdf5::save(model, filename, "tmp");

    std::vector<size_t> labeling = {1,1,1,1,1,1,1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0.0);

    ConstraintPool cp(200.0);
    cp.add_constraint(ConstraintPool::DetectionConstraint(2, 3));
    cp.add_constraint(ConstraintPool::OutgoingConstraint(3, 4, {5,6}));
    cp.add_constraint(ConstraintPool::IncomingConstraint({0,1}, 2));

    // save constraint pool to temp file
    std::string cp_filename(tmpnam(NULL));
    {
        std::ofstream cp_out_str(cp_filename.c_str());
        boost::archive::text_oarchive oa(cp_out_str);
        oa & cp;
    }

    opengm::ICM<OpengmModelDeprecated::ogmGraphicalModel, OpengmModelDeprecated::ogmAccumulator> inf(model);
    cp.add_constraints_to_problem(model, inf);

    // test wrong labelings
    labeling = {1, 1, 1, 1, 1, 1, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 200.0);

    labeling = {1, 2, 0, 1, 1, 2, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 400.0);

    // test good labelings
    labeling = {1, 1, 2, 2, 0, 0, 2};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);

    labeling = {0, 2, 2, 2, 0, 1, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);

    labeling = {1, 0, 1, 1, 1, 1, 1};
    BOOST_CHECK_EQUAL(model.evaluate(labeling.begin()), 0);

    //---------------------------------------------------------------------
    // load serialized model and constraint pool
    //---------------------------------------------------------------------
    OpengmModelDeprecated::ogmGraphicalModel model2;
    opengm::hdf5::load(model2, filename, "tmp");
    ConstraintPool cp2;
    {
        std::ifstream cp_in_str(cp_filename.c_str());
        boost::archive::text_iarchive ia(cp_in_str);
        ia & cp2;
    }

    labeling = {1,1,1,1,1,1,1};
    BOOST_CHECK_EQUAL(model2.evaluate(labeling.begin()), 0.0);

    opengm::ICM<OpengmModelDeprecated::ogmGraphicalModel, OpengmModelDeprecated::ogmAccumulator> inf2(model2);
    cp2.add_constraints_to_problem(model2, inf2);

    // test wrong labelings
    labeling = {1, 1, 1, 1, 1, 1, 1};
    BOOST_CHECK_EQUAL(model2.evaluate(labeling.begin()), 200.0);

    labeling = {1, 2, 0, 1, 1, 2, 1};
    BOOST_CHECK_EQUAL(model2.evaluate(labeling.begin()), 400.0);

    // test good labelings
    labeling = {1, 1, 2, 2, 0, 0, 2};
    BOOST_CHECK_EQUAL(model2.evaluate(labeling.begin()), 0);

    labeling = {0, 2, 2, 2, 0, 1, 1};
    BOOST_CHECK_EQUAL(model2.evaluate(labeling.begin()), 0);

    labeling = {1, 0, 1, 1, 1, 1, 1};
    BOOST_CHECK_EQUAL(model2.evaluate(labeling.begin()), 0);

    // delete tempfile
    remove(filename.c_str());
}

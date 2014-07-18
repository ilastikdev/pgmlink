#define BOOST_TEST_MODULE constraint_test

#include <vector>
#include <map>
#include <set>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/bind.hpp>

#include <opengm/inference/icm.hxx>

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
}

BOOST_AUTO_TEST_CASE(ConstraintPool_Size_Test)
{
    std::vector<size_t> dummy_vars = {1,2,3};

    ConstraintPool cp;
    cp.add_constraint(ConstraintPool::IncomingConstraint(dummy_vars.begin(), dummy_vars.end(), 4));
    cp.add_constraint(ConstraintPool::IncomingConstraint(dummy_vars.begin(), dummy_vars.end(), 4));
    cp.add_constraint(ConstraintPool::IncomingConstraint(dummy_vars.begin(), dummy_vars.end(), 4));

    BOOST_CHECK_EQUAL(cp.get_num_constraints(), 3);

    cp.add_constraint(ConstraintPool::OutgoingConstraint(1,2, dummy_vars.begin(), dummy_vars.end()));
    cp.add_constraint(ConstraintPool::OutgoingConstraint(1,2, dummy_vars.begin(), dummy_vars.end()));

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

    ConstraintPool cp;
    cp.set_big_m(200.0);
    std::vector<size_t> indices = {0,1};
    cp.add_constraint(ConstraintPool::IncomingConstraint(indices.begin(), indices.end(), 2));

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


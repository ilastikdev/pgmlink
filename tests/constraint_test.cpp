#define BOOST_TEST_MODULE constraint_test

#include <vector>
#include <map>
#include <set>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/bind.hpp>

#include "pgmlink/constraint_pool.hxx"

using namespace pgmlink::pgm;

BOOST_AUTO_TEST_CASE(IncomingFunctionTest)
{
    // T1, T2, V nodes with label space 3
    std::vector<size_t> shape = {3, 3, 3};

    IncomingConstraintFunction<double, size_t, size_t> constraint_func(shape.begin(), shape.end());
    constraint_func.set_forbidden_energy(200.0);
    std::vector<size_t> labeling = {1,2,3};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 0.0);

    labeling = {1,2,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 200.0);

    labeling = {1,0,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 0.0);

    labeling = {1,1,2};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 0.0);

    labeling = {1,1,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 200.0);
}

BOOST_AUTO_TEST_CASE(OutgoingFunctionTest)
{
    // A, D, T1, T2
    std::vector<size_t> shape = {3, 3, 3, 3};

    OutgoingConstraintFunction<double, size_t, size_t> constraint_func(shape.begin(), shape.end());
    constraint_func.set_forbidden_energy(200.0);
    std::vector<size_t> labeling = {3,0,2,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 0.0);

    labeling = {1,0,2,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 200.0);

    labeling = {1,1,1,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 0.0);

    labeling = {1,1,2,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 0.0);

    labeling = {2,1,1,2}; // division not allowed when A > 1!
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 200.0);
}

BOOST_AUTO_TEST_CASE(DetectionFunctionTest)
{
    // V,A nodes with label space 3
    std::vector<size_t> shape = {3, 3};

    DetectionConstraintFunction<double, size_t, size_t> constraint_func(shape.begin(), shape.end());
    constraint_func.set_forbidden_energy(200.0);
    std::vector<size_t> labeling = {1,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 0.0);

    labeling = {1,2};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 200.0);

    labeling = {1,0};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 0.0);

    labeling = {0,2};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 0.0);

    labeling = {2,1};
    BOOST_CHECK_EQUAL(constraint_func(labeling.begin(), labeling.end()), 200.0);
}

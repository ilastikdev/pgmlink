#ifndef CONSERVATIONTRACKING_PARAMETER_H
#define CONSERVATIONTRACKING_PARAMETER_H

#include <boost/function.hpp>
#include <boost/python.hpp>

#include "pgmlink/features/feature.h"
#include "pgmlink/traxels.h"
#include "pgmlink/uncertaintyParameter.h"

namespace pgmlink
{

enum class SolverType
{
    CplexSolver,
    DynProgSolver,
    DPInitCplexSolver
};

class Parameter
{
public:
    Parameter(
        unsigned int max_number_objects,
        boost::function<double (const Traxel&, const size_t)> detection,
        boost::function<double (const Traxel&, const size_t)> division,
        boost::function<double (const double)> transition,
        double forbidden_cost = 0,
        double ep_gap = 0.01,
        bool with_tracklets = false,
        bool with_divisions = true,
        boost::function<double (const Traxel&)> disappearance_cost_fn = ConstantFeature(500.0),
        boost::function<double (const Traxel&)> appearance_cost_fn = ConstantFeature(500.0),
        bool with_misdetections_allowed = true,
        bool with_appearance = true,
        bool with_disappearance = true,
        double transition_parameter = 5,
        bool with_constraints = true,
        UncertaintyParameter uncertainty_param = UncertaintyParameter(),
        double cplex_timeout = 1e75,
        double division_weight = 10,
        double detection_weight = 10,
        double transition_weight = 10,
        double border_width = 0,
        boost::python::object transition_classifier = boost::python::object(),
        bool with_optical_correction = false,
        SolverType solver = SolverType::CplexSolver,
        unsigned int num_threads = 0):
    max_number_objects(max_number_objects),
    detection(detection),
    division(division),
    transition(transition),
    border_width(border_width),
    forbidden_cost(forbidden_cost),
    ep_gap(ep_gap),
    with_tracklets(with_tracklets),
    with_divisions(with_divisions),
    disappearance_cost_fn(disappearance_cost_fn),
    appearance_cost_fn(appearance_cost_fn),
    with_misdetections_allowed(with_misdetections_allowed),
    with_appearance(with_appearance),
    with_disappearance(with_disappearance),
    transition_parameter(transition_parameter),
    with_constraints(with_constraints),
    uncertainty_param(uncertainty_param),
    cplex_timeout(cplex_timeout),
    division_weight(division_weight),
    detection_weight(detection_weight),
    transition_weight(transition_weight),
    transition_classifier(transition_classifier),
    with_optical_correction(with_optical_correction),
    solver(solver),
    num_threads(num_threads),
    with_swap(true),
    max_number_paths(std::numeric_limits<size_t>::max())
    {}

    // empty parameter needed for python
    Parameter() {}

    // settings
    unsigned int max_number_objects;
    boost::function<double (const Traxel&, const size_t)> detection;
    boost::function<double (const Traxel&, const size_t)> division;
    boost::function<double (const double)> transition;
    boost::function<double (const Traxel&, const Traxel&, const Traxel&)> motion_model3;
    boost::function<double (const Traxel&, const Traxel&, const Traxel&, const Traxel&)> motion_model4;
    double motion_model3_default;
    double motion_model4_default;
    double forbidden_cost;
    bool with_tracklets;
    bool with_divisions;
    boost::function<double (const Traxel&)> disappearance_cost_fn;
    boost::function<double (const Traxel&)> appearance_cost_fn;
    bool with_misdetections_allowed;
    bool with_appearance;
    bool with_disappearance;
    double transition_parameter;
    bool with_constraints;
    double division_weight;
    double detection_weight;
    double transition_weight;
    double border_width;
    bool with_optical_correction;
    SolverType solver;

    // cplex settings
    unsigned int num_threads;
    double ep_gap;
    double cplex_timeout;

    // dynprog settings
    size_t max_number_paths;
    bool with_swap;

    // perturbation settings
    UncertaintyParameter uncertainty_param;
    boost::python::object transition_classifier;

private:
    // python extensions:
    double python_caller_det_div(boost::python::object func, const Traxel& t, const size_t state)
    {
        assert(1 == PyCallable_Check(func.ptr()));
        // PyGILState_STATE pygilstate = PyGILState_Ensure();
        boost::python::object py_result = func(t, state);
        double result = boost::python::extract<double>(py_result);
        // PyGILState_Release(pygilstate);
        return result;
    }

    double python_caller_dis_appear(boost::python::object func, const Traxel& t)
    {
        assert(1 == PyCallable_Check(func.ptr()));
        // PyGILState_STATE pygilstate = PyGILState_Ensure();
        boost::python::object py_result = func(t);
        double result = boost::python::extract<double>(py_result);
        // PyGILState_Release(pygilstate);
        return result;
    }

    double python_caller_trans(boost::python::object func, double distance)
    {
        assert(1 == PyCallable_Check(func.ptr()));
        // PyGILState_STATE pygilstate = PyGILState_Ensure();
        boost::python::object py_result = func(distance);
        double result = boost::python::extract<double>(py_result);
        // PyGILState_Release(pygilstate);
        return result;
    }

    double python_caller_motion_model3(boost::python::object func, const Traxel& a, const Traxel& b, const Traxel& c)
    {
        assert(1 == PyCallable_Check(func.ptr()));
        boost::python::object py_result = func(a, b, c);
        double result = boost::python::extract<double>(py_result);
        return result;
    }

    double python_caller_motion_model4(boost::python::object func, const Traxel& a, const Traxel& b, const Traxel& c, const Traxel& d)
    {
        assert(1 == PyCallable_Check(func.ptr()));
        boost::python::object py_result = func(a, b, c, d);
        double result = boost::python::extract<double>(py_result);
        return result;
    }
public:
    /// Expects a function with signature (Traxel traxel, size_t state) -> double energy
    void register_detection_func(boost::python::object func)
    {
        detection = boost::bind(&Parameter::python_caller_det_div, this, func, _1, _2);
    }

    /// Expects a function with signature (Traxel traxel, size_t state) -> double energy
    void register_division_func(boost::python::object func)
    {
        division = boost::bind(&Parameter::python_caller_det_div, this, func, _1, _2);
    }

    /// Expects a function with signature (double distance) -> double energy
    void register_transition_func(boost::python::object func)
    {
        transition = boost::bind(&Parameter::python_caller_trans, this, func, _1);
    }

    /// Expects a function with signature (Traxel traxel) -> double energy
    void register_appearance_func(boost::python::object func)
    {
        appearance_cost_fn = boost::bind(&Parameter::python_caller_dis_appear, this, func, _1);
    }

    /// Expects a function with signature (Traxel traxel) -> double energy
    void register_disappearance_func(boost::python::object func)
    {
        disappearance_cost_fn = boost::bind(&Parameter::python_caller_dis_appear, this, func, _1);
    }

    /// Expects a function with signature (Traxel, Traxel, Traxel) -> double energy
    void register_motion_model3_func(boost::python::object func, double default_value)
    {
        motion_model3 = boost::bind(&Parameter::python_caller_motion_model3, this, func, _1, _2, _3);
        motion_model3_default = default_value;
    }

    /// Expects a function with signature (Traxel, Traxel, Traxel, Traxel) -> double energy
    void register_motion_model4_func(boost::python::object func, double default_value)
    {
        motion_model4 = boost::bind(&Parameter::python_caller_motion_model4, this, func, _1, _2, _3, _4);
        motion_model4_default = default_value;
    }
};

} // end namespace pgmlink

#endif

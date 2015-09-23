#ifndef CONSTRAINT_POOL_EXPLICIT_HXX
#define CONSTRAINT_POOL_EXPLICIT_HXX

#include <map>
#include <vector>
#include <memory.h>
#include <string.h>
#include <opengm/inference/lpcplex.hxx>
#include <boost/serialization/serialization.hpp>

#include "pgmlink/pgm.h"
#include "pgmlink/constraint_function.hxx"

namespace pgmlink
{
namespace pgm
{

//typedef OpengmModelDeprecated::ogmGraphicalModel ConstraintPoolExplicitOpengmModel;
//typedef opengm::LPCplex<ConstraintPoolExplicitOpengmModel,OpengmModelDeprecated::ogmAccumulator> ConstraintPoolExplicitCplexOptimizer;
typedef PertExplicitGmType ConstraintPoolExplicitOpengmModel;
typedef opengm::LPCplex<pgmlink::PertExplicitGmType, OpengmModelDeprecated::ogmAccumulator> ConstraintPoolExplicitCplexOptimizer;

//------------------------------------------------------------------------
// ConstraintPoolExplicit
//------------------------------------------------------------------------
/// Of the incoming and outgoing constraint type, we will have instatiations of different order.
/// The constraint pool contains a map indexed by the "constraint-function-signatures",
/// pointing to the function object. When adding functions to the model we look up
/// whether it is already present in here and just get the function identifier to create a factor.
class ConstraintPoolExplicit
{
public:
    typedef ConstraintPoolExplicitOpengmModel::ValueType ValueType;
    typedef ConstraintPoolExplicitOpengmModel::LabelType LabelType;
    typedef ConstraintPoolExplicitOpengmModel::IndexType IndexType;

public:
    ConstraintPoolExplicit(ValueType big_m = 200.0,
                   bool with_divisions = true,
                   bool with_appearance = true,
                   bool with_disappearance = true,
                   bool with_misdetections = true):
        big_m_(big_m),
        with_divisions_(with_divisions),
        with_appearance_(with_appearance),
        with_disappearance_(with_disappearance),
        with_misdetections_(with_misdetections)
    {}

    template<class CONSTRAINT_TYPE>
    void add_constraint(const CONSTRAINT_TYPE& constraint);

    /// This method either adds all the factors to the graphical model
    /// or - given INF is opengm::LpCplex and "not force_softconstraint" - adds hard constraints
    /// to the CPLEX optimizer.
    template<class GM, class INF>
    void add_constraints_to_problem(GM& model, INF& optimizer);

    /// This method adds all the factors to the graphical model
    template<class GM, class INF>
    void add_constraints_to_model(GM& model, INF& optimizer);

    /// Allow to add this set of constraints to a different model with the given index mapping (key is index of variable in original model)
    /// Only those constraints or functions are instanciated where all participating indices do have a mapping
    template<class GM, class INF>
    void add_constraints_to_problem(GM& model, INF& optimizer, std::map<size_t, size_t>& index_mapping);

    template<class GM, class INF>
    void add_constraints_to_model(GM& model, INF& optimizer, std::map<size_t, size_t>& index_mapping);

    size_t get_num_constraints()
    {
        return incoming_constraints_.size() + outgoing_constraints_.size() + detection_constraints_.size();
    }

    void force_softconstraint(bool enable)
    {
        if(enable)
        {
            LOG(logWARNING) << "[ConstraintPoolExplicit]: Forcing soft constraint not yet implemented";
        }
        force_softconstraint_ = enable;
    }

    void set_big_m(ValueType m)
    {
        big_m_ = m;
    }

public:
    // types to store constraint instanciations
    class IncomingConstraint
    {
    public:
        IncomingConstraint(const std::vector<IndexType>& transition_nodes,
                           IndexType disappearance_node):
            transition_nodes(transition_nodes),
            disappearance_node(disappearance_node)
        {}

        std::vector<IndexType> transition_nodes;
        IndexType disappearance_node;

    private:
        // boost serialization interface
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int);
        IncomingConstraint() {}
    };

    class IncomingLinearConstraint
    {
    public:
        IncomingLinearConstraint(const std::vector<IndexType>& transition_nodes,
                           IndexType disappearance_node):
            transition_nodes(transition_nodes),
            disappearance_node(disappearance_node)
        {}

        std::vector<IndexType> transition_nodes;
        IndexType disappearance_node;

    private:
        // boost serialization interface
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int);
        IncomingLinearConstraint() {}
    };

    class OutgoingConstraint
    {
    public:
        OutgoingConstraint(IndexType appearance_node,
                           int division_node,
                           const std::vector<IndexType>& transition_nodes):
            appearance_node(appearance_node),
            division_node(division_node),
            transition_nodes(transition_nodes)
        {}

        IndexType appearance_node;
        int division_node;
        std::vector<IndexType> transition_nodes;

    private:
        // boost serialization interface
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int);
        OutgoingConstraint() {}
    };

    class OutgoingLinearConstraint
    {
    public:
        OutgoingLinearConstraint(IndexType appearance_node,
                           int division_node,
                           const std::vector<IndexType>& transition_nodes):
            appearance_node(appearance_node),
            division_node(division_node),
            transition_nodes(transition_nodes)
        {}

        IndexType appearance_node;
        int division_node;
        std::vector<IndexType> transition_nodes;

    private:
        // boost serialization interface
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int);
        OutgoingLinearConstraint() {}
    };

    class DetectionConstraint
    {
    public:
        DetectionConstraint(IndexType disappearance_node,
                            IndexType appearance_node):
            disappearance_node(disappearance_node),
            appearance_node(appearance_node)
        {}

        IndexType disappearance_node;
        IndexType appearance_node;

    private:
        // boost serialization interface
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int);
        DetectionConstraint() {}
    };

    class DetectionLinearConstraint
    {
    public:
        DetectionLinearConstraint(IndexType disappearance_node,
                            IndexType appearance_node):
            disappearance_node(disappearance_node),
            appearance_node(appearance_node)
        {}

        IndexType disappearance_node;
        IndexType appearance_node;

    private:
        // boost serialization interface
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int);
        DetectionLinearConstraint() {}
    };

    class FixNodeValueConstraint
    {
    public:
        FixNodeValueConstraint(IndexType node,
                               size_t value):
            node(node),
            value(value)
        {}

        IndexType node;
        size_t value;

    private:
        // boost serialization interface
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int);
        FixNodeValueConstraint() {}
    };

    class FixNodeValueLinearConstraint
    {
    public:
        FixNodeValueLinearConstraint(IndexType node,
                               size_t value):
            node(node),
            value(value)
        {}

        IndexType node;
        size_t value;

    private:
        // boost serialization interface
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int);
        FixNodeValueLinearConstraint() {}
    };

protected:
    template<class CONSTRAINT_TYPE>
    void constraint_indices(std::vector<IndexType>& indices, const CONSTRAINT_TYPE& constraint);

    template<class GM, class INF, class FUNCTION_TYPE, class CONSTRAINT_TYPE>
    void add_constraint_type_to_problem(GM& model, INF&, const std::vector<CONSTRAINT_TYPE>& constraints);

    template<class GM, class INF, class FUNCTION_TYPE, class CONSTRAINT_TYPE>
    void add_constraint_type_to_model(GM& model, INF&, const std::vector<CONSTRAINT_TYPE>& constraints);

    template<class FUNCTION_TYPE, class CONSTRAINT_TYPE>
    void configure_function(FUNCTION_TYPE* func, CONSTRAINT_TYPE constraint);

    template<class CONSTRAINT_TYPE>
    bool check_all_constraint_vars_in_mapping(std::map<size_t, size_t>& index_mapping, const CONSTRAINT_TYPE& constraint);
protected:
    std::vector<IncomingConstraint> incoming_constraints_;
    std::vector<OutgoingConstraint> outgoing_constraints_;
    std::vector<OutgoingConstraint> outgoing_no_div_constraints_;
    std::vector<DetectionConstraint> detection_constraints_;
    std::vector<FixNodeValueConstraint> fix_node_value_constraints_;

    std::vector<IncomingLinearConstraint> incoming_linear_constraints_;
    std::vector<OutgoingLinearConstraint> outgoing_linear_constraints_;
    std::vector<OutgoingLinearConstraint> outgoing_no_div_linear_constraints_;
    std::vector<DetectionLinearConstraint> detection_linear_constraints_;
    std::vector<FixNodeValueLinearConstraint> fix_node_value_linear_constraints_;

    ValueType big_m_;
    bool with_divisions_;
    bool with_misdetections_;
    bool with_appearance_;
    bool with_disappearance_;

    bool force_softconstraint_;

private:
    // boost serialization interface
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int);
};


//------------------------------------------------------------------------
// ConstraintPoolExplicit - Implementation
//------------------------------------------------------------------------

template<class CONSTRAINT_TYPE>
void ConstraintPoolExplicit::add_constraint(const CONSTRAINT_TYPE&)
{
    throw std::logic_error("only template specializations of this method should be called");
}

template<>
void ConstraintPoolExplicit::add_constraint(const ConstraintPoolExplicit::IncomingConstraint& constraint);

template<>
void ConstraintPoolExplicit::add_constraint(const ConstraintPoolExplicit::OutgoingConstraint& constraint);

template<>
void ConstraintPoolExplicit::add_constraint(const ConstraintPoolExplicit::DetectionConstraint& constraint);

template<>
void ConstraintPoolExplicit::add_constraint(const ConstraintPoolExplicit::FixNodeValueConstraint& constraint);

template<>
void ConstraintPoolExplicit::add_constraint(const ConstraintPoolExplicit::IncomingLinearConstraint& constraint);

template<>
void ConstraintPoolExplicit::add_constraint(const ConstraintPoolExplicit::OutgoingLinearConstraint& constraint);

template<>
void ConstraintPoolExplicit::add_constraint(const ConstraintPoolExplicit::DetectionLinearConstraint& constraint);

template<>
void ConstraintPoolExplicit::add_constraint(const ConstraintPoolExplicit::FixNodeValueLinearConstraint& constraint);

//------------------------------------------------------------------------
template<class GM, class INF>
void ConstraintPoolExplicit::add_constraints_to_problem(GM& model, INF& inf)
{
    // TODO: handle force_softconstraint_

    std::cout << "===========================================================>in explicit add_constraints_to_problem" << std::endl;
    add_constraint_type_to_problem<GM, INF, IncomingConstraintFunction<ValueType, IndexType, LabelType>, IncomingConstraint>(model, inf, incoming_constraints_);
    add_constraint_type_to_problem<GM, INF, OutgoingConstraintFunction<ValueType, IndexType, LabelType>, OutgoingConstraint>(model, inf, outgoing_constraints_);
    add_constraint_type_to_problem<GM, INF, OutgoingNoDivConstraintFunction<ValueType, IndexType, LabelType>, OutgoingConstraint>(model, inf, outgoing_no_div_constraints_);
    add_constraint_type_to_problem<GM, INF, DetectionConstraintFunction<ValueType, IndexType, LabelType>, DetectionConstraint>(model, inf, detection_constraints_);
    add_constraint_type_to_problem<GM, INF, FixNodeValueConstraintFunction<ValueType, IndexType, LabelType>, FixNodeValueConstraint>(model, inf, fix_node_value_constraints_);
}

template<class GM, class INF>
void ConstraintPoolExplicit::add_constraints_to_model(GM& model, INF& inf)
{
    // TODO: handle force_softconstraint_

    std::cout << "=========@@@==================================================>in explicit add_constraints_to_model" << std::endl;
    add_constraint_type_to_model<GM, INF, IncomingLinearConstraintFunction<ValueType, IndexType, LabelType>, IncomingLinearConstraint>(model, inf, incoming_linear_constraints_);
    add_constraint_type_to_model<GM, INF, OutgoingLinearConstraintFunction<ValueType, IndexType, LabelType>, OutgoingLinearConstraint>(model, inf, outgoing_linear_constraints_);
    add_constraint_type_to_model<GM, INF, OutgoingNoDivLinearConstraintFunction<ValueType, IndexType, LabelType>, OutgoingLinearConstraint>(model, inf, outgoing_no_div_linear_constraints_);
    add_constraint_type_to_model<GM, INF, DetectionLinearConstraintFunction<ValueType, IndexType, LabelType>, DetectionLinearConstraint>(model, inf, detection_linear_constraints_);
    add_constraint_type_to_model<GM, INF, FixNodeValueLinearConstraintFunction<ValueType, IndexType, LabelType>, FixNodeValueLinearConstraint>(model, inf, fix_node_value_linear_constraints_);
    std::cout << "=========@@@==================================================>out explicit add_constraints_to_model" << std::endl;
}

template<class GM, class INF>
void ConstraintPoolExplicit::add_constraints_to_problem(GM& model, INF& inf, std::map<size_t, size_t>& index_mapping)
{
    std::cout << "===========================================================>in explicit add_constraints_to_problem" << std::endl;
    std::vector<IncomingConstraint> remapped_incoming_constraints;
    std::vector<OutgoingConstraint> remapped_outgoing_constraints;
    std::vector<OutgoingConstraint> remapped_outgoing_no_div_constraints;
    std::vector<DetectionConstraint> remapped_detection_constraints;
    std::vector<FixNodeValueConstraint> remapped_fix_node_value_constraints;

    for(auto constraint : incoming_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
        {
            continue;
        }

        size_t disappearance_node = index_mapping[constraint.disappearance_node];
        std::vector<size_t> transition_nodes;
        for(auto t : constraint.transition_nodes)
        {
            transition_nodes.push_back(index_mapping[t]);
        }
        remapped_incoming_constraints.push_back(IncomingConstraint(transition_nodes, disappearance_node));
    }

    for(auto constraint : outgoing_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
        {
            continue;
        }

        size_t appearance_node = index_mapping[constraint.appearance_node];
        size_t division_node = index_mapping[constraint.division_node];
        std::vector<size_t> transition_nodes;
        for(auto t : constraint.transition_nodes)
        {
            transition_nodes.push_back(index_mapping[t]);
        }
        remapped_outgoing_constraints.push_back(OutgoingConstraint(appearance_node, division_node, transition_nodes));
    }

    for(auto constraint : outgoing_no_div_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
        {
            continue;
        }

        size_t appearance_node = index_mapping[constraint.appearance_node];
        std::vector<size_t> transition_nodes;
        for(auto t : constraint.transition_nodes)
        {
            transition_nodes.push_back(index_mapping[t]);
        }
        remapped_outgoing_no_div_constraints.push_back(OutgoingConstraint(appearance_node, -1, transition_nodes));
    }

    for(auto constraint : detection_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
        {
            continue;
        }

        size_t appearance_node = index_mapping[constraint.appearance_node];
        size_t disappearance_node = index_mapping[constraint.disappearance_node];
        remapped_detection_constraints.push_back(DetectionConstraint(disappearance_node, appearance_node));
    }

    for(auto constraint : fix_node_value_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
        {
            continue;
        }

        size_t node = index_mapping[constraint.node];
        remapped_fix_node_value_constraints.push_back(FixNodeValueConstraint(node, constraint.value));
    }

    add_constraint_type_to_problem<GM, INF, IncomingConstraintFunction<ValueType, IndexType, LabelType>, IncomingConstraint>(model, inf, remapped_incoming_constraints);
    add_constraint_type_to_problem<GM, INF, OutgoingConstraintFunction<ValueType, IndexType, LabelType>, OutgoingConstraint>(model, inf, remapped_outgoing_constraints);
    add_constraint_type_to_problem<GM, INF, OutgoingNoDivConstraintFunction<ValueType, IndexType, LabelType>, OutgoingConstraint>(model, inf, remapped_outgoing_no_div_constraints);
    add_constraint_type_to_problem<GM, INF, DetectionConstraintFunction<ValueType, IndexType, LabelType>, DetectionConstraint>(model, inf, remapped_detection_constraints);
    add_constraint_type_to_problem<GM, INF, FixNodeValueConstraintFunction<ValueType, IndexType, LabelType>, FixNodeValueConstraint>(model, inf, remapped_fix_node_value_constraints);
}

template<class GM, class INF>
void ConstraintPoolExplicit::add_constraints_to_model(GM& model, INF& inf, std::map<size_t, size_t>& index_mapping)
{
    std::cout << "=====0======================================================>in explicit add_constraints_to_model" << std::endl;
    std::vector<IncomingLinearConstraint> remapped_incoming_linear_constraints;
    std::vector<OutgoingLinearConstraint> remapped_outgoing_linear_constraints;
    std::vector<OutgoingLinearConstraint> remapped_outgoing_no_div_linear_constraints;
    std::vector<DetectionLinearConstraint> remapped_detection_linear_constraints;
    std::vector<FixNodeValueLinearConstraint> remapped_fix_node_value_linear_constraints;

    for(auto constraint : incoming_linear_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
        {
            continue;
        }

        size_t disappearance_node = index_mapping[constraint.disappearance_node];
        std::vector<size_t> transition_nodes;
        for(auto t : constraint.transition_nodes)
        {
            transition_nodes.push_back(index_mapping[t]);
        }
        remapped_incoming_linear_constraints.push_back(IncomingLinearConstraint(transition_nodes, disappearance_node));
    }

    for(auto constraint : outgoing_linear_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
        {
            continue;
        }

        size_t appearance_node = index_mapping[constraint.appearance_node];
        size_t division_node = index_mapping[constraint.division_node];
        std::vector<size_t> transition_nodes;
        for(auto t : constraint.transition_nodes)
        {
            transition_nodes.push_back(index_mapping[t]);
        }
        remapped_outgoing_linear_constraints.push_back(OutgoingLinearConstraint(appearance_node, division_node, transition_nodes));
    }

    for(auto constraint : outgoing_no_div_linear_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
        {
            continue;
        }

        size_t appearance_node = index_mapping[constraint.appearance_node];
        std::vector<size_t> transition_nodes;
        for(auto t : constraint.transition_nodes)
        {
            transition_nodes.push_back(index_mapping[t]);
        }
        remapped_outgoing_no_div_linear_constraints.push_back(OutgoingLinearConstraint(appearance_node, -1, transition_nodes));
    }

    for(auto constraint : detection_linear_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
        {
            continue;
        }

        size_t appearance_node = index_mapping[constraint.appearance_node];
        size_t disappearance_node = index_mapping[constraint.disappearance_node];
        remapped_detection_linear_constraints.push_back(DetectionLinearConstraint(disappearance_node, appearance_node));
    }
/*
    for(auto constraint : fix_node_value_linear_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
        {
            continue;
        }

        size_t node = index_mapping[constraint.node];
        remapped_fix_node_value_linear_constraints.push_back(FixNodeValueLinearConstraint(node, constraint.value));
    }
    */
    std::cout << "=====1======================================================>in explicit add_constraints_to_model" << std::endl;

    add_constraint_type_to_model<GM, INF, IncomingLinearConstraintFunction<ValueType, IndexType, LabelType>, IncomingLinearConstraint>(model, inf, remapped_incoming_linear_constraints);
    add_constraint_type_to_model<GM, INF, OutgoingLinearConstraintFunction<ValueType, IndexType, LabelType>, OutgoingLinearConstraint>(model, inf, remapped_outgoing_linear_constraints);
    add_constraint_type_to_model<GM, INF, OutgoingNoDivLinearConstraintFunction<ValueType, IndexType, LabelType>, OutgoingLinearConstraint>(model, inf, remapped_outgoing_no_div_linear_constraints);
    add_constraint_type_to_model<GM, INF, DetectionLinearConstraintFunction<ValueType, IndexType, LabelType>, DetectionLinearConstraint>(model, inf, remapped_detection_linear_constraints);
    //add_constraint_type_to_model<GM, INF, FixNodeValueLinearConstraintFunction<ValueType, IndexType, LabelType>, FixNodeValueLinearConstraint>(model, inf, remapped_fix_node_value_linear_constraints);
    std::cout << "=====2======================================================>out explicit add_constraints_to_model" << std::endl;
}

template<class GM, class INF, class FUNCTION_TYPE, class CONSTRAINT_TYPE>
void ConstraintPoolExplicit::add_constraint_type_to_problem(GM& model, INF&, const std::vector<CONSTRAINT_TYPE>& constraints)
{
    std::cout << "===========================================================>in explicit add_constraint_type_to_problem" << std::endl;
    LOG(logINFO) << "[ConstraintPoolExplicit]: Using soft constraints";
    std::map< std::vector<IndexType>, FUNCTION_TYPE* > constraint_functions;
    for(typename std::vector<CONSTRAINT_TYPE>::const_iterator it = constraints.begin(); it != constraints.end(); ++it)
    {
        const CONSTRAINT_TYPE& constraint = *it;
        std::vector<IndexType> indices;
        constraint_indices(indices, constraint);
        if(indices.size() < 2)
        {
            continue;
        }

        std::vector<IndexType> shape;
        for(std::vector<IndexType>::iterator idx = indices.begin(); idx != indices.end(); ++idx)
        {
            shape.push_back(model.numberOfLabels(*idx));
        }

        // see if the function is already present in our map and model
        if(constraint_functions.find(shape) == constraint_functions.end())
        {
            constraint_functions[shape] = new FUNCTION_TYPE(shape.begin(), shape.end());
            constraint_functions[shape]->set_forbidden_energy(big_m_);
            configure_function(constraint_functions[shape], *it);
        }

        // create factor
        OpengmFactor< FUNCTION_TYPE > factor(*constraint_functions[shape], indices.begin(), indices.end());

        std::cout << "===========================================================>ADDING FACTOR TO MODEL" << std::endl;
        factor.add_to(model);
    }
}

template<class GM, class INF, class FUNCTION_TYPE, class CONSTRAINT_TYPE>
void ConstraintPoolExplicit::add_constraint_type_to_model(GM& model, INF&, const std::vector<CONSTRAINT_TYPE>& constraints)
{
    std::cout << "===========================================================>in explicit add_constraint_type_to_model" << std::endl;
    LOG(logINFO) << "[ConstraintPoolExplicit]: Using soft constraints";
    std::map< std::vector<IndexType>, FUNCTION_TYPE* > constraint_functions;
    for(typename std::vector<CONSTRAINT_TYPE>::const_iterator it = constraints.begin(); it != constraints.end(); ++it)
    {
        const CONSTRAINT_TYPE& constraint = *it;
        std::vector<IndexType> indices;
        constraint_indices(indices, constraint);
        if(indices.size() < 2)
        {
            continue;
        }

        std::vector<IndexType> shape;
        for(std::vector<IndexType>::iterator idx = indices.begin(); idx != indices.end(); ++idx)
        {
            shape.push_back(model.numberOfLabels(*idx));
        }

        // see if the function is already present in our map and model
        if(constraint_functions.find(shape) == constraint_functions.end())
        {
            constraint_functions[shape] = new FUNCTION_TYPE(shape.begin(), shape.end());
            constraint_functions[shape]->set_forbidden_energy(big_m_);
            configure_function(constraint_functions[shape], *it);
        }

        // create factor
        OpengmFactor< FUNCTION_TYPE > factor(*constraint_functions[shape], indices.begin(), indices.end());

        std::cout << "=====-------======================================================>ADDING CONSTRAINT TO MODEL" << std::endl;
        factor.add_to(model);
    }
}

//------------------------------------------------------------------------
// specialization for IncomingConstraintFunction
template<>
void ConstraintPoolExplicit::add_constraint_type_to_problem<ConstraintPoolExplicitOpengmModel,
     ConstraintPoolExplicitCplexOptimizer,
     IncomingConstraintFunction<ConstraintPoolExplicit::ValueType, ConstraintPoolExplicit::IndexType, ConstraintPoolExplicit::LabelType>,
     ConstraintPoolExplicit::IncomingConstraint>
     (
         ConstraintPoolExplicitOpengmModel& model,
         ConstraintPoolExplicitCplexOptimizer& optimizer,
         const std::vector<ConstraintPoolExplicit::IncomingConstraint>& constraints
     );

//------------------------------------------------------------------------
// specialization for OutgoingConstraintFunction
template<>
void ConstraintPoolExplicit::add_constraint_type_to_problem<ConstraintPoolExplicitOpengmModel,
     ConstraintPoolExplicitCplexOptimizer,
     OutgoingConstraintFunction<ConstraintPoolExplicit::ValueType, ConstraintPoolExplicit::IndexType, ConstraintPoolExplicit::LabelType>,
     ConstraintPoolExplicit::OutgoingConstraint>
     (
         ConstraintPoolExplicitOpengmModel& model,
         ConstraintPoolExplicitCplexOptimizer& optimizer,
         const std::vector<ConstraintPoolExplicit::OutgoingConstraint>& constraints
     );

template<>
void ConstraintPoolExplicit::add_constraint_type_to_problem<ConstraintPoolExplicitOpengmModel,
     ConstraintPoolExplicitCplexOptimizer,
     OutgoingNoDivConstraintFunction<ConstraintPoolExplicit::ValueType, ConstraintPoolExplicit::IndexType, ConstraintPoolExplicit::LabelType>,
     ConstraintPoolExplicit::OutgoingConstraint>
     (
         ConstraintPoolExplicitOpengmModel& model,
         ConstraintPoolExplicitCplexOptimizer& optimizer,
         const std::vector<ConstraintPoolExplicit::OutgoingConstraint>& constraints
     );

//------------------------------------------------------------------------
// specialization for DetectionConstraintFunction
template<>
void ConstraintPoolExplicit::add_constraint_type_to_problem<ConstraintPoolExplicitOpengmModel,
     ConstraintPoolExplicitCplexOptimizer,
     DetectionConstraintFunction<ConstraintPoolExplicit::ValueType, ConstraintPoolExplicit::IndexType, ConstraintPoolExplicit::LabelType>,
     ConstraintPoolExplicit::DetectionConstraint>
     (
         ConstraintPoolExplicitOpengmModel& model,
         ConstraintPoolExplicitCplexOptimizer& optimizer,
         const std::vector<ConstraintPoolExplicit::DetectionConstraint>& constraints
     );

//------------------------------------------------------------------------
// specialization for FixNodeValueConstraintFunction
template<>
void ConstraintPoolExplicit::add_constraint_type_to_problem<ConstraintPoolExplicitOpengmModel,
     ConstraintPoolExplicitCplexOptimizer,
     FixNodeValueConstraintFunction<ConstraintPoolExplicit::ValueType, ConstraintPoolExplicit::IndexType, ConstraintPoolExplicit::LabelType>,
     ConstraintPoolExplicit::FixNodeValueConstraint>
     (
         ConstraintPoolExplicitOpengmModel& model,
         ConstraintPoolExplicitCplexOptimizer& optimizer,
         const std::vector<ConstraintPoolExplicit::FixNodeValueConstraint>& constraints
     );

//------------------------------------------------------------------------
// specialization for IncomingLinearConstraintFunction
template<>
void ConstraintPoolExplicit::add_constraint_type_to_model<ConstraintPoolExplicitOpengmModel,
     ConstraintPoolExplicitCplexOptimizer,
     IncomingLinearConstraintFunction<ConstraintPoolExplicit::ValueType, ConstraintPoolExplicit::IndexType, ConstraintPoolExplicit::LabelType>,
     ConstraintPoolExplicit::IncomingLinearConstraint>
     (
         ConstraintPoolExplicitOpengmModel& model,
         ConstraintPoolExplicitCplexOptimizer& optimizer,
         const std::vector<ConstraintPoolExplicit::IncomingLinearConstraint>& constraints
     );

//------------------------------------------------------------------------
// specialization for OutgoingLinearConstraintFunction
template<>
void ConstraintPoolExplicit::add_constraint_type_to_model<ConstraintPoolExplicitOpengmModel,
     ConstraintPoolExplicitCplexOptimizer,
     OutgoingLinearConstraintFunction<ConstraintPoolExplicit::ValueType, ConstraintPoolExplicit::IndexType, ConstraintPoolExplicit::LabelType>,
     ConstraintPoolExplicit::OutgoingLinearConstraint>
     (
         ConstraintPoolExplicitOpengmModel& model,
         ConstraintPoolExplicitCplexOptimizer& optimizer,
         const std::vector<ConstraintPoolExplicit::OutgoingLinearConstraint>& constraints
     );

template<>
void ConstraintPoolExplicit::add_constraint_type_to_model<ConstraintPoolExplicitOpengmModel,
     ConstraintPoolExplicitCplexOptimizer,
     OutgoingNoDivLinearConstraintFunction<ConstraintPoolExplicit::ValueType, ConstraintPoolExplicit::IndexType, ConstraintPoolExplicit::LabelType>,
     ConstraintPoolExplicit::OutgoingLinearConstraint>
     (
         ConstraintPoolExplicitOpengmModel& model,
         ConstraintPoolExplicitCplexOptimizer& optimizer,
         const std::vector<ConstraintPoolExplicit::OutgoingLinearConstraint>& constraints
     );

//------------------------------------------------------------------------
// specialization for DetectionLinearConstraintFunction
template<>
void ConstraintPoolExplicit::add_constraint_type_to_model<ConstraintPoolExplicitOpengmModel,
     ConstraintPoolExplicitCplexOptimizer,
     DetectionLinearConstraintFunction<ConstraintPoolExplicit::ValueType, ConstraintPoolExplicit::IndexType, ConstraintPoolExplicit::LabelType>,
     ConstraintPoolExplicit::DetectionLinearConstraint>
     (
         ConstraintPoolExplicitOpengmModel& model,
         ConstraintPoolExplicitCplexOptimizer& optimizer,
         const std::vector<ConstraintPoolExplicit::DetectionLinearConstraint>& constraints
     );

//------------------------------------------------------------------------
// specialization for FixNodeValueLinearConstraintFunction
template<>
void ConstraintPoolExplicit::add_constraint_type_to_model<ConstraintPoolExplicitOpengmModel,
     ConstraintPoolExplicitCplexOptimizer,
     FixNodeValueLinearConstraintFunction<ConstraintPoolExplicit::ValueType, ConstraintPoolExplicit::IndexType, ConstraintPoolExplicit::LabelType>,
     ConstraintPoolExplicit::FixNodeValueLinearConstraint>
     (
         ConstraintPoolExplicitOpengmModel& model,
         ConstraintPoolExplicitCplexOptimizer& optimizer,
         const std::vector<ConstraintPoolExplicit::FixNodeValueLinearConstraint>& constraints
     );

//------------------------------------------------------------------------
template<class CONSTRAINT_TYPE>
void ConstraintPoolExplicit::constraint_indices(std::vector<ConstraintPoolExplicit::IndexType>&, const CONSTRAINT_TYPE&)
{
    throw std::logic_error("only template specializations of this method should be called");
}

template<>
void ConstraintPoolExplicit::constraint_indices<ConstraintPoolExplicit::IncomingConstraint>(std::vector<ConstraintPoolExplicit::IndexType>& indices, const IncomingConstraint& constraint);

template<>
void ConstraintPoolExplicit::constraint_indices(std::vector<ConstraintPoolExplicit::IndexType>& indices, const OutgoingConstraint& constraint);

template<>
void ConstraintPoolExplicit::constraint_indices(std::vector<ConstraintPoolExplicit::IndexType>& indices, const DetectionConstraint& constraint);

template<>
void ConstraintPoolExplicit::constraint_indices(std::vector<ConstraintPoolExplicit::IndexType>& indices, const FixNodeValueConstraint& constraint);

template<>
void ConstraintPoolExplicit::constraint_indices<ConstraintPoolExplicit::IncomingLinearConstraint>(std::vector<ConstraintPoolExplicit::IndexType>& indices, const IncomingLinearConstraint& constraint);

template<>
void ConstraintPoolExplicit::constraint_indices(std::vector<ConstraintPoolExplicit::IndexType>& indices, const OutgoingLinearConstraint& constraint);

template<>
void ConstraintPoolExplicit::constraint_indices(std::vector<ConstraintPoolExplicit::IndexType>& indices, const DetectionLinearConstraint& constraint);

template<>
void ConstraintPoolExplicit::constraint_indices(std::vector<ConstraintPoolExplicit::IndexType>& indices, const FixNodeValueLinearConstraint& constraint);

//------------------------------------------------------------------------
template<class CONSTRAINT_TYPE>
bool ConstraintPoolExplicit::check_all_constraint_vars_in_mapping(std::map<size_t, size_t>& index_mapping, const CONSTRAINT_TYPE& constraint)
{
    std::vector<size_t> indices;
    constraint_indices(indices, constraint);

    for(size_t i = 0; i < indices.size(); i++)
    {
        if(index_mapping.find(indices[i]) == index_mapping.end())
        {
            return false;
        }
    }

    return true;
}

//------------------------------------------------------------------------
template<class FUNCTION_TYPE, class CONSTRAINT_TYPE>
void ConstraintPoolExplicit::configure_function(FUNCTION_TYPE* func, CONSTRAINT_TYPE constraint)
{
    throw std::logic_error("only template specializations of this method should be called");
}

template<>
void ConstraintPoolExplicit::configure_function(IncomingConstraintFunction<ValueType, IndexType, LabelType>*, ConstraintPoolExplicit::IncomingConstraint);

template<>
void ConstraintPoolExplicit::configure_function(OutgoingNoDivConstraintFunction<ValueType, IndexType, LabelType>*, ConstraintPoolExplicit::OutgoingConstraint);

template<>
void ConstraintPoolExplicit::configure_function(OutgoingConstraintFunction<ValueType, IndexType, LabelType>* func, ConstraintPoolExplicit::IncomingConstraint);

template<>
void ConstraintPoolExplicit::configure_function(DetectionConstraintFunction<ValueType, IndexType, LabelType>* func, ConstraintPoolExplicit::DetectionConstraint);

template<>
void ConstraintPoolExplicit::configure_function(FixNodeValueConstraintFunction<ValueType, IndexType, LabelType>* func, ConstraintPoolExplicit::FixNodeValueConstraint constraint);

template<>
void ConstraintPoolExplicit::configure_function(IncomingLinearConstraintFunction<ValueType, IndexType, LabelType>*, ConstraintPoolExplicit::IncomingLinearConstraint);

template<>
void ConstraintPoolExplicit::configure_function(OutgoingNoDivLinearConstraintFunction<ValueType, IndexType, LabelType>*, ConstraintPoolExplicit::OutgoingLinearConstraint);

template<>
void ConstraintPoolExplicit::configure_function(OutgoingLinearConstraintFunction<ValueType, IndexType, LabelType>* func, ConstraintPoolExplicit::IncomingLinearConstraint);

template<>
void ConstraintPoolExplicit::configure_function(DetectionLinearConstraintFunction<ValueType, IndexType, LabelType>* func, ConstraintPoolExplicit::DetectionLinearConstraint);

template<>
void ConstraintPoolExplicit::configure_function(FixNodeValueLinearConstraintFunction<ValueType, IndexType, LabelType>* func, ConstraintPoolExplicit::FixNodeValueLinearConstraint constraint);

//------------------------------------------------------------------------
// Serialization
//------------------------------------------------------------------------
template<class Archive>
void ConstraintPoolExplicit::serialize(Archive & ar, const unsigned int)
{
    // config first
    ar & big_m_;
    ar & with_divisions_;
    ar & with_misdetections_;
    ar & with_appearance_;
    ar & with_disappearance_;
    ar & force_softconstraint_;

    // constraint arrays
    ar & incoming_constraints_;
    ar & outgoing_constraints_;
    ar & outgoing_no_div_constraints_;
    ar & detection_constraints_;
    ar & fix_node_value_constraints_;
}

template<class Archive>
void ConstraintPoolExplicit::IncomingConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & transition_nodes;
    ar & disappearance_node;
}

template<class Archive>
void ConstraintPoolExplicit::OutgoingConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & appearance_node;
    ar & division_node;
    ar & transition_nodes;
}

template<class Archive>
void ConstraintPoolExplicit::DetectionConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & appearance_node;
    ar & disappearance_node;
}

template<class Archive>
void ConstraintPoolExplicit::FixNodeValueConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & value;
    ar & node;
}

template<class Archive>
void ConstraintPoolExplicit::IncomingLinearConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & transition_nodes;
    ar & disappearance_node;
}

template<class Archive>
void ConstraintPoolExplicit::OutgoingLinearConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & appearance_node;
    ar & division_node;
    ar & transition_nodes;
}

template<class Archive>
void ConstraintPoolExplicit::DetectionLinearConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & appearance_node;
    ar & disappearance_node;
}

template<class Archive>
void ConstraintPoolExplicit::FixNodeValueLinearConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & value;
    ar & node;
}

} // namespace pgm
} // namespace pgmlink

#endif // CONSTRAINT_POOL_EXPLICIT_HXX

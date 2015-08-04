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

    /// Allow to add this set of constraints to a different model with the given index mapping (key is index of variable in original model)
    /// Only those constraints or functions are instanciated where all participating indices do have a mapping
    template<class GM, class INF>
    void add_constraints_to_problem(GM& model, INF& optimizer, std::map<size_t, size_t>& index_mapping);

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

protected:
    template<class CONSTRAINT_TYPE>
    void constraint_indices(std::vector<IndexType>& indices, const CONSTRAINT_TYPE& constraint);

    template<class GM, class INF, class FUNCTION_TYPE, class CONSTRAINT_TYPE>
    void add_constraint_type_to_problem(GM& model, INF&, const std::vector<CONSTRAINT_TYPE>& constraints);

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

//------------------------------------------------------------------------
template<class GM, class INF>
void ConstraintPoolExplicit::add_constraints_to_problem(GM& model, INF& inf)
{
    // TODO: handle force_softconstraint_

    add_constraint_type_to_problem<GM, INF, IncomingConstraintFunction<ValueType, IndexType, LabelType>, IncomingConstraint>(model, inf, incoming_constraints_);
    add_constraint_type_to_problem<GM, INF, OutgoingConstraintFunction<ValueType, IndexType, LabelType>, OutgoingConstraint>(model, inf, outgoing_constraints_);
    add_constraint_type_to_problem<GM, INF, OutgoingNoDivConstraintFunction<ValueType, IndexType, LabelType>, OutgoingConstraint>(model, inf, outgoing_no_div_constraints_);
    add_constraint_type_to_problem<GM, INF, DetectionConstraintFunction<ValueType, IndexType, LabelType>, DetectionConstraint>(model, inf, detection_constraints_);
    add_constraint_type_to_problem<GM, INF, FixNodeValueConstraintFunction<ValueType, IndexType, LabelType>, FixNodeValueConstraint>(model, inf, fix_node_value_constraints_);
}

template<class GM, class INF>
void ConstraintPoolExplicit::add_constraints_to_problem(GM& model, INF& inf, std::map<size_t, size_t>& index_mapping)
{
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
    add_constraint_type_to_problem<GM, INF, FixNodeValueConstraintFunction<ValueType, IndexType, LabelType>, FixNodeValueConstraint>(model, inf, remapped_detection_constraints);
}

template<class GM, class INF, class FUNCTION_TYPE, class CONSTRAINT_TYPE>
void ConstraintPoolExplicit::add_constraint_type_to_problem(GM& model, INF&, const std::vector<CONSTRAINT_TYPE>& constraints)
{
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

} // namespace pgm
} // namespace pgmlink

#endif // CONSTRAINT_POOL_EXPLICIT_HXX

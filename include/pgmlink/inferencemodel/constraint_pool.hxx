#ifndef CONSTRAINT_POOL_HXX
#define CONSTRAINT_POOL_HXX

#include <map>
#include <vector>
#include <memory.h>
#include <string.h>
#include <opengm/inference/lpcplex.hxx>
#include <boost/serialization/serialization.hpp>

#include "pgm.h"
#include "constraint_function.hxx"

namespace pgmlink
{
namespace pgm
{

//typedef OpengmModelDeprecated::ogmGraphicalModel ConstraintPoolOpengmModel;
//typedef opengm::LPCplex<ConstraintPoolConstraintPoolOpengmModel,OpengmModelDeprecated::ogmAccumulator> ConstraintPoolCplexOptimizer;
typedef PertGmType ConstraintPoolOpengmModel;
typedef opengm::LPCplex<pgmlink::PertGmType,OpengmModelDeprecated::ogmAccumulator> ConstraintPoolCplexOptimizer;

//------------------------------------------------------------------------
// ConstraintPool
//------------------------------------------------------------------------
/// Of the incoming and outgoing constraint type, we will have instatiations of different order.
/// The constraint pool contains a map indexed by the "constraint-function-signatures",
/// pointing to the function object. When adding functions to the model we look up
/// whether it is already present in here and just get the function identifier to create a factor.
class ConstraintPool
{
public:
    typedef ConstraintPoolOpengmModel::ValueType ValueType;
    typedef ConstraintPoolOpengmModel::LabelType LabelType;
    typedef ConstraintPoolOpengmModel::IndexType IndexType;

public:
    ConstraintPool(ValueType big_m = 200.0,
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
            LOG(logWARNING) << "[ConstraintPool]: Forcing soft constraint not yet implemented";
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
        IncomingConstraint(){}
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

protected:
    template<class CONSTRAINT_TYPE>
    void constraint_indices(std::vector<IndexType>& indices, const CONSTRAINT_TYPE& constraint);

    template<class GM, class INF, class FUNCTION_TYPE, class CONSTRAINT_TYPE>
    void add_constraint_type_to_problem(GM& model, INF&, const std::vector<CONSTRAINT_TYPE>& constraints);

    template<class FUNCTION_TYPE>
    void configure_function(FUNCTION_TYPE* func);

    template<class CONSTRAINT_TYPE>
    bool check_all_constraint_vars_in_mapping(std::map<size_t, size_t>& index_mapping, const CONSTRAINT_TYPE& constraint);
protected:
    std::vector<IncomingConstraint> incoming_constraints_;
    std::vector<OutgoingConstraint> outgoing_constraints_;
    std::vector<OutgoingConstraint> outgoing_no_div_constraints_;
    std::vector<DetectionConstraint> detection_constraints_;

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
// ConstraintPool - Implementation
//------------------------------------------------------------------------

template<class CONSTRAINT_TYPE>
void ConstraintPool::add_constraint(const CONSTRAINT_TYPE&)
{
    throw std::logic_error("only template specializations of this method should be called");
}

template<>
void ConstraintPool::add_constraint(const ConstraintPool::IncomingConstraint& constraint);

template<>
void ConstraintPool::add_constraint(const ConstraintPool::OutgoingConstraint& constraint);

template<>
void ConstraintPool::add_constraint(const ConstraintPool::DetectionConstraint& constraint);

//------------------------------------------------------------------------
template<class GM, class INF>
void ConstraintPool::add_constraints_to_problem(GM& model, INF& inf)
{
    // TODO: handle force_softconstraint_

    add_constraint_type_to_problem<GM, INF, IncomingConstraintFunction<ValueType,IndexType,LabelType>, IncomingConstraint>(model, inf, incoming_constraints_);
    add_constraint_type_to_problem<GM, INF, OutgoingConstraintFunction<ValueType,IndexType,LabelType>, OutgoingConstraint>(model, inf, outgoing_constraints_);
    add_constraint_type_to_problem<GM, INF, OutgoingNoDivConstraintFunction<ValueType,IndexType,LabelType>, OutgoingConstraint>(model, inf, outgoing_no_div_constraints_);
    add_constraint_type_to_problem<GM, INF, DetectionConstraintFunction<ValueType,IndexType,LabelType>, DetectionConstraint>(model, inf, detection_constraints_);
}

template<class GM, class INF>
void ConstraintPool::add_constraints_to_problem(GM& model, INF& inf, std::map<size_t, size_t>& index_mapping)
{
    std::vector<IncomingConstraint> remapped_incoming_constraints;
    std::vector<OutgoingConstraint> remapped_outgoing_constraints;
    std::vector<OutgoingConstraint> remapped_outgoing_no_div_constraints;
    std::vector<DetectionConstraint> remapped_detection_constraints;

    for(auto constraint : incoming_constraints_)
    {
        if(!check_all_constraint_vars_in_mapping(index_mapping, constraint))
            continue;

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
            continue;

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
            continue;

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
            continue;

        size_t appearance_node = index_mapping[constraint.appearance_node];
        size_t disappearance_node = index_mapping[constraint.disappearance_node];
        remapped_detection_constraints.push_back(DetectionConstraint(disappearance_node, appearance_node));
    }

    add_constraint_type_to_problem<GM, INF, IncomingConstraintFunction<ValueType,IndexType,LabelType>, IncomingConstraint>(model, inf, remapped_incoming_constraints);
    add_constraint_type_to_problem<GM, INF, OutgoingConstraintFunction<ValueType,IndexType,LabelType>, OutgoingConstraint>(model, inf, remapped_outgoing_constraints);
    add_constraint_type_to_problem<GM, INF, OutgoingNoDivConstraintFunction<ValueType,IndexType,LabelType>, OutgoingConstraint>(model, inf, remapped_outgoing_no_div_constraints);
    add_constraint_type_to_problem<GM, INF, DetectionConstraintFunction<ValueType,IndexType,LabelType>, DetectionConstraint>(model, inf, remapped_detection_constraints);
}

template<class GM, class INF, class FUNCTION_TYPE, class CONSTRAINT_TYPE>
void ConstraintPool::add_constraint_type_to_problem(GM& model, INF&, const std::vector<CONSTRAINT_TYPE>& constraints)
{
    LOG(logINFO) << "[ConstraintPool]: Using soft constraints";
    std::map< std::vector<IndexType>, FUNCTION_TYPE* > constraint_functions;
    for(typename std::vector<CONSTRAINT_TYPE>::const_iterator it = constraints.begin(); it != constraints.end(); ++it)
    {
        const CONSTRAINT_TYPE& constraint = *it;
        std::vector<IndexType> indices;
        constraint_indices(indices, constraint);
        if(indices.size() < 2)
            continue;

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
            configure_function(constraint_functions[shape]);
        }

        // create factor
        OpengmFactor< FUNCTION_TYPE > factor(*constraint_functions[shape], indices.begin(), indices.end());
        factor.add_to(model);
    }
}

//------------------------------------------------------------------------
// specialization for IncomingConstraintFunction
template<>
void ConstraintPool::add_constraint_type_to_problem<ConstraintPoolOpengmModel,
ConstraintPoolCplexOptimizer,
IncomingConstraintFunction<ConstraintPool::ValueType, ConstraintPool::IndexType, ConstraintPool::LabelType>,
ConstraintPool::IncomingConstraint>
(
        ConstraintPoolOpengmModel& model,
        ConstraintPoolCplexOptimizer& optimizer,
        const std::vector<ConstraintPool::IncomingConstraint>& constraints
);

//------------------------------------------------------------------------
// specialization for OutgoingConstraintFunction
template<>
void ConstraintPool::add_constraint_type_to_problem<ConstraintPoolOpengmModel,
ConstraintPoolCplexOptimizer,
OutgoingConstraintFunction<ConstraintPool::ValueType, ConstraintPool::IndexType, ConstraintPool::LabelType>,
ConstraintPool::OutgoingConstraint>
(
        ConstraintPoolOpengmModel& model,
        ConstraintPoolCplexOptimizer& optimizer,
        const std::vector<ConstraintPool::OutgoingConstraint>& constraints
);

template<>
void ConstraintPool::add_constraint_type_to_problem<ConstraintPoolOpengmModel,
ConstraintPoolCplexOptimizer,
OutgoingNoDivConstraintFunction<ConstraintPool::ValueType, ConstraintPool::IndexType, ConstraintPool::LabelType>,
ConstraintPool::OutgoingConstraint>
(
        ConstraintPoolOpengmModel& model,
        ConstraintPoolCplexOptimizer& optimizer,
        const std::vector<ConstraintPool::OutgoingConstraint>& constraints
);

//------------------------------------------------------------------------
// specialization for DetectionConstraintFunction
template<>
void ConstraintPool::add_constraint_type_to_problem<ConstraintPoolOpengmModel,
ConstraintPoolCplexOptimizer,
DetectionConstraintFunction<ConstraintPool::ValueType, ConstraintPool::IndexType, ConstraintPool::LabelType>,
ConstraintPool::DetectionConstraint>
(
        ConstraintPoolOpengmModel& model,
        ConstraintPoolCplexOptimizer& optimizer,
        const std::vector<ConstraintPool::DetectionConstraint>& constraints
);

//------------------------------------------------------------------------
template<class CONSTRAINT_TYPE>
void ConstraintPool::constraint_indices(std::vector<ConstraintPool::IndexType>&, const CONSTRAINT_TYPE&)
{
    throw std::logic_error("only template specializations of this method should be called");
}

template<>
void ConstraintPool::constraint_indices<ConstraintPool::IncomingConstraint>(std::vector<ConstraintPool::IndexType>& indices, const IncomingConstraint& constraint);

template<>
void ConstraintPool::constraint_indices(std::vector<ConstraintPool::IndexType>& indices, const OutgoingConstraint& constraint);

template<>
void ConstraintPool::constraint_indices(std::vector<ConstraintPool::IndexType>& indices, const DetectionConstraint& constraint);

//------------------------------------------------------------------------
template<class CONSTRAINT_TYPE>
bool ConstraintPool::check_all_constraint_vars_in_mapping(std::map<size_t, size_t>& index_mapping, const CONSTRAINT_TYPE& constraint)
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
template<class FUNCTION_TYPE>
void ConstraintPool::configure_function(FUNCTION_TYPE* func)
{
    throw std::logic_error("only template specializations of this method should be called");
}

template<>
void ConstraintPool::configure_function(IncomingConstraintFunction<ValueType, IndexType, LabelType>*);

template<>
void ConstraintPool::configure_function(OutgoingNoDivConstraintFunction<ValueType, IndexType, LabelType>*);

template<>
void ConstraintPool::configure_function(OutgoingConstraintFunction<ValueType, IndexType, LabelType>* func);

template<>
void ConstraintPool::configure_function(DetectionConstraintFunction<ValueType, IndexType, LabelType>* func);

//------------------------------------------------------------------------
// Serialization
//------------------------------------------------------------------------
template<class Archive>
void ConstraintPool::serialize(Archive & ar, const unsigned int)
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
}

template<class Archive>
void ConstraintPool::IncomingConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & transition_nodes;
    ar & disappearance_node;
}

template<class Archive>
void ConstraintPool::OutgoingConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & appearance_node;
    ar & division_node;
    ar & transition_nodes;
}


template<class Archive>
void ConstraintPool::DetectionConstraint::serialize(Archive & ar, const unsigned int)
{
    ar & appearance_node;
    ar & disappearance_node;
}

} // namespace pgm
} // namespace pgmlink

#endif // CONSTRAINT_POOL_HXX

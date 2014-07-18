#ifndef CONSTRAINT_POOL_HXX
#define CONSTRAINT_POOL_HXX

#include <map>
#include <vector>
#include <opengm/inference/lpcplex.hxx>

#include "pgm.h"
#include "constraint_function.hxx"

namespace pgmlink
{
namespace pgm
{

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
    typedef OpengmModelDeprecated::ogmGraphicalModel::ValueType ValueType;
    typedef OpengmModelDeprecated::ogmGraphicalModel::LabelType LabelType;
    typedef OpengmModelDeprecated::ogmGraphicalModel::IndexType IndexType;

public:
    template<class CONSTRAINT_TYPE>
    void add_constraint(const CONSTRAINT_TYPE& constraint);

    /// This method either adds all the factors to the graphical model
    /// or - given INF is opengm::LpCplex and "not force_hardconstraint" - adds hard constraints
    /// to the CPLEX optimizer.
    template<class GM, class INF>
    void add_constraints_to_problem(GM& model,
                                    INF& optimizer);

    void set_big_m(ValueType big_m)
    {
        big_m_ = big_m;
    }

    size_t get_num_constraints()
    {
        return incoming_constraints_.size() + outgoing_constraints_.size() + detection_constraints_.size();
    }

public:
    // types to store constraint instanciations
    struct IncomingConstraint
    {
        template<class TRANS_VAR_ITERATOR>
        IncomingConstraint(TRANS_VAR_ITERATOR transition_nodes_begin,
                           TRANS_VAR_ITERATOR transition_nodes_end,
                           IndexType disappearance_node):
            transition_nodes(transition_nodes_begin, transition_nodes_end),
            disappearance_node(disappearance_node)
        {}

        std::vector<IndexType> transition_nodes;
        IndexType disappearance_node;
    };

    struct OutgoingConstraint
    {
        template<class TRANS_VAR_ITERATOR>
        OutgoingConstraint(IndexType appearance_node,
                           IndexType division_node,
                           TRANS_VAR_ITERATOR transition_nodes_begin,
                           TRANS_VAR_ITERATOR transition_nodes_end):
            appearance_node(appearance_node),
            division_node(division_node),
            transition_nodes(transition_nodes_begin, transition_nodes_end)
        {}

        IndexType appearance_node;
        IndexType division_node;
        std::vector<IndexType> transition_nodes;
    };

    struct DetectionConstraint
    {
        DetectionConstraint(IndexType disappearance_node,
                            IndexType appearance_node):
            disappearance_node(disappearance_node),
            appearance_node(appearance_node)
        {}

        IndexType disappearance_node;
        IndexType appearance_node;
    };

protected:
    void instanciate_incoming_constraints();
    void instanciate_outgoing_constraints();
    void instanciate_detection_constraints();

    template<class CONSTRAINT_TYPE>
    void constraint_indices(std::vector<IndexType>& indices, const CONSTRAINT_TYPE& constraint);

    template<class GM, class INF, class FUNCTION_TYPE, class CONSTRAINT_TYPE>
    void add_constraint_type_to_problem(GM& model, INF&, const std::vector<CONSTRAINT_TYPE>& constraints);

protected:
    ValueType big_m_;
    std::vector<IncomingConstraint> incoming_constraints_;
    std::vector<OutgoingConstraint> outgoing_constraints_;
    std::vector<DetectionConstraint> detection_constraints_;
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
void ConstraintPool::add_constraint(const ConstraintPool::IncomingConstraint& constraint)
{
    incoming_constraints_.push_back(constraint);
}

template<>
void ConstraintPool::add_constraint(const ConstraintPool::OutgoingConstraint& constraint)
{
    outgoing_constraints_.push_back(constraint);
}

template<>
void ConstraintPool::add_constraint(const ConstraintPool::DetectionConstraint& constraint)
{
    detection_constraints_.push_back(constraint);
}

//------------------------------------------------------------------------
template<class GM, class INF>
void ConstraintPool::add_constraints_to_problem(GM& model, INF& inf)
{
    add_constraint_type_to_problem<GM, INF, IncomingConstraintFunction<ValueType,IndexType,LabelType>, IncomingConstraint>(model, inf, incoming_constraints_);
    add_constraint_type_to_problem<GM, INF, OutgoingConstraintFunction<ValueType,IndexType,LabelType>, OutgoingConstraint>(model, inf, outgoing_constraints_);
    add_constraint_type_to_problem<GM, INF, DetectionConstraintFunction<ValueType,IndexType,LabelType>, DetectionConstraint>(model, inf, detection_constraints_);
}

template<class GM, class INF, class FUNCTION_TYPE, class CONSTRAINT_TYPE>
void ConstraintPool::add_constraint_type_to_problem(GM& model, INF&, const std::vector<CONSTRAINT_TYPE>& constraints)
{
    std::map< std::vector<IndexType>, FUNCTION_TYPE* > constraint_functions;
    for(typename std::vector<CONSTRAINT_TYPE>::const_iterator it = constraints.begin(); it != constraints.end(); ++it)
    {
        const CONSTRAINT_TYPE& constraint = *it;
        std::vector<IndexType> indices;
        constraint_indices(indices, constraint);

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
        }

        // create factor
        OpengmFactor< FUNCTION_TYPE > factor(*constraint_functions[shape], indices.begin(), indices.end());
        factor.add_to(model);
    }
}

//------------------------------------------------------------------------
template<class CONSTRAINT_TYPE>
void ConstraintPool::constraint_indices(std::vector<ConstraintPool::IndexType>& indices, const CONSTRAINT_TYPE& constraint)
{
    throw std::logic_error("only template specializations of this method should be called");
}

template<>
void ConstraintPool::constraint_indices<ConstraintPool::IncomingConstraint>(std::vector<ConstraintPool::IndexType>& indices, const IncomingConstraint& constraint)
{
    indices.insert(indices.begin(), constraint.transition_nodes.begin(), constraint.transition_nodes.end());
    indices.push_back(constraint.disappearance_node);
}

template<>
void ConstraintPool::constraint_indices(std::vector<ConstraintPool::IndexType>& indices, const OutgoingConstraint& constraint)
{
    indices.push_back(constraint.appearance_node);
    indices.push_back(constraint.division_node);
    indices.insert(indices.begin(), constraint.transition_nodes.begin(), constraint.transition_nodes.end());
}

template<>
void ConstraintPool::constraint_indices(std::vector<ConstraintPool::IndexType>& indices, const DetectionConstraint& constraint)
{
    indices.push_back(constraint.disappearance_node);
    indices.push_back(constraint.appearance_node);
}

//// member function template specialization to handle CPLEX
//template<>
//void ConstraintPool::add_constraints_to_problem(OpengmModelDeprecated::ogmGraphicalModel& model,
//                                                       opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator>& optimizer)
//{
//    std::cout << "Trying to instanciate constraints" << std::endl;
//    throw std::logic_error("not yet implemented");
//}

// TODO: add member function template specialization for Gurobi

} // namespace pgm
} // namespace pgmlink

#endif // CONSTRAINT_POOL_HXX

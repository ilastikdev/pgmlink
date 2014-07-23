#ifndef CONSTRAINT_POOL_HXX
#define CONSTRAINT_POOL_HXX

#include <map>
#include <vector>
#include <opengm/inference/lpcplex.hxx>
#include <boost/serialization/serialization.hpp>

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
    void add_constraints_to_problem(GM& model,
                                    INF& optimizer);

    size_t get_num_constraints()
    {
        return incoming_constraints_.size() + outgoing_constraints_.size() + detection_constraints_.size();
    }

    void force_softconstraint(bool enable)
    {
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
void ConstraintPool::add_constraint(const ConstraintPool::IncomingConstraint& constraint)
{
    incoming_constraints_.push_back(constraint);
}

template<>
void ConstraintPool::add_constraint(const ConstraintPool::OutgoingConstraint& constraint)
{
    // here we separate the outgoing constraints with division node from those without,
    // such that the template specializations work
    if(constraint.division_node >= 0)
        outgoing_constraints_.push_back(constraint);
    else
        outgoing_no_div_constraints_.push_back(constraint);
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
    // TODO: handle force_softconstraint_

    add_constraint_type_to_problem<GM, INF, IncomingConstraintFunction<ValueType,IndexType,LabelType>, IncomingConstraint>(model, inf, incoming_constraints_);
    add_constraint_type_to_problem<GM, INF, OutgoingConstraintFunction<ValueType,IndexType,LabelType>, OutgoingConstraint>(model, inf, outgoing_constraints_);
    add_constraint_type_to_problem<GM, INF, OutgoingNoDivConstraintFunction<ValueType,IndexType,LabelType>, OutgoingConstraint>(model, inf, outgoing_no_div_constraints_);
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
void ConstraintPool::add_constraint_type_to_problem<OpengmModelDeprecated::ogmGraphicalModel,
opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator>,
IncomingConstraintFunction<ConstraintPool::ValueType, ConstraintPool::IndexType, ConstraintPool::LabelType>,
ConstraintPool::IncomingConstraint>
(
        OpengmModelDeprecated::ogmGraphicalModel& model,
        opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator>& optimizer,
        const std::vector<ConstraintPool::IncomingConstraint>& constraints
)
{
    for(auto it = constraints.begin(); it != constraints.end(); ++it)
    {
        const ConstraintPool::IncomingConstraint& constraint = *it;

        // nothing to do if no incoming
        if(constraint.transition_nodes.size() == 0)
            continue;

        // 0 <= sum_nu [ nu * sum_i (Y_ij[nu] ) ] - sum_nu ( nu * X_j[nu] ) - sum_nu ( nu * Dis_j[nu] ) <= 0
        std::vector<size_t> cplex_idxs;
        std::vector<int> coeffs;
        std::stringstream constraint_name;
        constraint_name << "incoming transitions: sum( transition-nodes ";

        for (auto incoming_it = constraint.transition_nodes.begin(); incoming_it != constraint.transition_nodes.end(); ++incoming_it)
        {
            for (size_t state = 1; state < model.numberOfLabels(*incoming_it); ++state)
            {
                cplex_idxs.push_back(optimizer.lpNodeVi(*incoming_it, state));
                coeffs.push_back(state);
            }
            constraint_name << *incoming_it << " ";
        }

        constraint_name << ") = disappearance-node " << constraint.disappearance_node;

        for (size_t state = 1; state < model.numberOfLabels(constraint.disappearance_node); ++state)
        {
            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.disappearance_node, state));
            coeffs.push_back(-state);
        }

        optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 0, constraint_name.str().c_str());
        LOG(logDEBUG3) << constraint_name.str();
    }
}

//------------------------------------------------------------------------
// specialization for OutgoingConstraintFunction
template<>
void ConstraintPool::add_constraint_type_to_problem<OpengmModelDeprecated::ogmGraphicalModel,
opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator>,
OutgoingConstraintFunction<ConstraintPool::ValueType, ConstraintPool::IndexType, ConstraintPool::LabelType>,
ConstraintPool::OutgoingConstraint>
(
        OpengmModelDeprecated::ogmGraphicalModel& model,
        opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator>& optimizer,
        const std::vector<ConstraintPool::OutgoingConstraint>& constraints
)
{
    for(auto it = constraints.begin(); it != constraints.end(); ++it)
    {
        const ConstraintPool::OutgoingConstraint& constraint = *it;

        // nothing to do if no outgoing
        if(constraint.transition_nodes.size() == 0)
            continue;

        std::vector<size_t> cplex_idxs, cplex_idxs2;
        std::vector<int> coeffs, coeffs2;

        std::stringstream constraint_name;
        constraint_name << "outgoing transitions: sum( transition-nodes ";

        // couple detection and transitions: Y_ij <= App_i
        for (auto outgoing_it = constraint.transition_nodes.begin(); outgoing_it != constraint.transition_nodes.end(); ++outgoing_it)
        {
            for(size_t t_state = 0; t_state < model.numberOfLabels(*outgoing_it) - 1; t_state++)
            {
                for(size_t a_state = t_state + 1; a_state < model.numberOfLabels(constraint.appearance_node); a_state++)
                {
                    cplex_idxs.clear();
                    coeffs.clear();
                    coeffs.push_back(1);
                    cplex_idxs2.push_back(optimizer.lpNodeVi(constraint.appearance_node, a_state));
                    coeffs.push_back(1);
                    cplex_idxs.push_back(optimizer.lpNodeVi(*outgoing_it, t_state));
                    constraint_name.str(std::string()); // clear the name
                    constraint_name << "outgoing: 0 <= App_i[" << a_state << "] + Y_ij[" << t_state << "] <= 1; ";
                    constraint_name << "g.id(n) = " << *outgoing_it << ", g.id(a) = " << constraint.appearance_node;
                    optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(),
                            0, 1, constraint_name.str().c_str());
                    LOG(logDEBUG3) << constraint_name.str();
                }
            }
        }


        int div_cplex_id = -1;
        if (with_divisions_ && constraint.division_node >= 0)
        {
            LOG(logDEBUG3) << "div_node_map_[n] = " << constraint.division_node;
            LOG(logDEBUG3) << "number_of_transition_nodes_ = " << constraint.transition_nodes.size();
            // LOG(logDEBUG3) << "number_of_division_nodes_ = " << number_of_division_nodes_; ???
            div_cplex_id = optimizer.lpNodeVi(constraint.division_node, 1);
        }

        if (constraint.transition_nodes.size() > 0)
        {
            // couple transitions: sum(Y_ij) = D_i + App_i
            cplex_idxs.clear();
            coeffs.clear();

            for (auto outgoing_it = constraint.transition_nodes.begin(); outgoing_it != constraint.transition_nodes.end(); ++outgoing_it)
            {
                for (size_t state = 1; state < model.numberOfLabels(*outgoing_it); ++state)
                {
                    coeffs.push_back(state);
                    cplex_idxs.push_back(optimizer.lpNodeVi(*outgoing_it, state));
                }
            }

            if (div_cplex_id != -1)
            {
                cplex_idxs.push_back(div_cplex_id);
                coeffs.push_back(-1);
            }

            for (size_t state = 1; state < model.numberOfLabels(constraint.appearance_node); ++state)
            {
                coeffs.push_back(-state);
                cplex_idxs.push_back(optimizer.lpNodeVi(constraint.appearance_node, state));
            }

            // 0 <= sum_nu [ sum_j( nu * Y_ij[nu] ) ] - [ sum_nu nu * X_i[nu] + D_i[1] + sum_nu nu * App_i[nu] ]<= 0
            constraint_name.str(std::string()); // clear the name
            constraint_name << "couple transitions: ";
            constraint_name << " sum(Y_ij) = D_i + App_i added for nodes: App_i=" << constraint.appearance_node
                            << ", D_i = " << constraint.division_node;
            optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 0,
                    constraint_name.str().c_str());
            LOG(logDEBUG3) << constraint_name.str();
        }

        if (div_cplex_id != -1)
        {
            // couple detection and division: D_i = 1 => App_i = 1
            cplex_idxs.clear();
            coeffs.clear();

            cplex_idxs.push_back(div_cplex_id);
            coeffs.push_back(1);

            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.appearance_node, 1));
            coeffs.push_back(-1);

            // -1 <= D_i[1] - App_i[1] <= 0
            constraint_name.str(std::string()); // clear the name
            constraint_name << "couple division and detection: ";
            constraint_name << " D_i=1 => App_i =1 added for n = "
                    << constraint.appearance_node << ", d = " << constraint.division_node;
            optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), -1, 0,
                    constraint_name.str().c_str());
            LOG(logDEBUG3) << constraint_name.str();

            // couple divsion and transition: D_1 = 1 => sum_k(Y_ik) = 2
            cplex_idxs2.clear();
            coeffs2.clear(); // -m <= 2 * D_i[1] - sum_j ( Y_ij[1] ) <= 0
            cplex_idxs2.push_back(div_cplex_id);
            coeffs2.push_back(2);

            for (auto outgoing_it = constraint.transition_nodes.begin(); outgoing_it != constraint.transition_nodes.end(); ++outgoing_it)
            {
                for (size_t state = 2; state < model.numberOfLabels(*outgoing_it); ++state)
                {
                    // D_i[1] = 1 => Y_ij[nu] = 0 forall nu > 1
                    cplex_idxs.clear();
                    coeffs.clear();
                    cplex_idxs.push_back(div_cplex_id);
                    coeffs.push_back(1);
                    cplex_idxs.push_back(optimizer.lpNodeVi(*outgoing_it, state));
                    coeffs.push_back(1);

                    // 0 <= D_i[1] + Y_ij[nu] <= 1 forall nu>1
                    constraint_name.str(std::string()); // clear the name
                    constraint_name << "couple division and transition: ";
                    constraint_name << " D_i=1 => Y_i[nu]=0 added for "
                            << "d = " << constraint.division_node << ", y = " << *outgoing_it << ", nu = "
                            << state;

                    optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(),
                            0, 1, constraint_name.str().c_str());
                    LOG(logDEBUG3) << constraint_name.str();

                }

                cplex_idxs2.push_back(optimizer.lpNodeVi(*outgoing_it, 1));
                coeffs2.push_back(-1);
            }

            // -m <= 2 * D_i[1] - sum_j (Y_ij[1]) <= 0
            constraint_name.str(std::string()); // clear the name
            constraint_name << "couple division and transitions: ";
            constraint_name  << " D_i = 1 => sum_k(Y_ik) = 2 added for "
                    << "d = " << constraint.division_node;
            optimizer.addConstraint(cplex_idxs2.begin(), cplex_idxs2.end(), coeffs2.begin(),
                    -int(model.numberOfLabels(constraint.appearance_node)-1), 0, constraint_name.str().c_str());
            LOG(logDEBUG3) << constraint_name.str();
        }
    }
}

template<>
void ConstraintPool::add_constraint_type_to_problem<OpengmModelDeprecated::ogmGraphicalModel,
opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator>,
OutgoingNoDivConstraintFunction<ConstraintPool::ValueType, ConstraintPool::IndexType, ConstraintPool::LabelType>,
ConstraintPool::OutgoingConstraint>
(
        OpengmModelDeprecated::ogmGraphicalModel& model,
        opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator>& optimizer,
        const std::vector<ConstraintPool::OutgoingConstraint>& constraints
)
{
    // for the CPLEX specialization we do the same for with and without divisions
    add_constraint_type_to_problem<OpengmModelDeprecated::ogmGraphicalModel,
                                    opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator>,
                                    OutgoingConstraintFunction<ValueType,IndexType,LabelType>,
                                    OutgoingConstraint>
            (model, optimizer, constraints);
}

//------------------------------------------------------------------------
// specialization for DetectionConstraintFunction
template<>
void ConstraintPool::add_constraint_type_to_problem<OpengmModelDeprecated::ogmGraphicalModel,
opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator>,
DetectionConstraintFunction<ConstraintPool::ValueType, ConstraintPool::IndexType, ConstraintPool::LabelType>,
ConstraintPool::DetectionConstraint>
(
        OpengmModelDeprecated::ogmGraphicalModel& model,
        opengm::LPCplex<OpengmModelDeprecated::ogmGraphicalModel,OpengmModelDeprecated::ogmAccumulator>& optimizer,
        const std::vector<ConstraintPool::DetectionConstraint>& constraints
)
{
    for(auto it = constraints.begin(); it != constraints.end(); ++it)
    {
        const ConstraintPool::DetectionConstraint& constraint = *it;

        // 0 <= sum_nu [ nu * sum_i (Y_ij[nu] ) ] - sum_nu ( nu * X_j[nu] ) - sum_nu ( nu * Dis_j[nu] ) <= 0
        std::vector<size_t> cplex_idxs;
        std::vector<int> coeffs;
        std::stringstream constraint_name;

        for (size_t state = 1; state < model.numberOfLabels(constraint.appearance_node); ++state)
        {
            cplex_idxs.clear();
            coeffs.clear();

            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.appearance_node, state));
            coeffs.push_back(1);

            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.disappearance_node, state));
            coeffs.push_back(-1);

            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.disappearance_node, 0));
            coeffs.push_back(-1);

            // A_i[nu] = 1 => V_i[nu] = 1 v V_i[0] = 1
            // -1 <= App_i[nu] - ( Dis_i[nu] + Dis_i[0] ) <= 0 forall nu > 0
            constraint_name.str(std::string()); // clear the name
            constraint_name << "disappearance/appearance coupling: ";
            constraint_name << " A_i[nu] = 1 => V_i[nu] = 1 v V_i[0] = 1 added for nodes "
                    << constraint.appearance_node << ", " << constraint.disappearance_node;
            constraint_name << " for state: " << state;
            optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), -1,
                    0, constraint_name.str().c_str());
            LOG(logDEBUG3) << constraint_name.str();
        }

        for (size_t state = 1; state < model.numberOfLabels(constraint.disappearance_node); ++state)
        {
            cplex_idxs.clear();
            coeffs.clear();

            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.disappearance_node, state));
            coeffs.push_back(1);

            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.appearance_node, state));
            coeffs.push_back(-1);

            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.appearance_node, 0));
            coeffs.push_back(-1);

            // V_i[nu] = 1 => A_i[nu] = 1 v A_i[0] = 1
            // -1 <= Dis_i[nu] - ( App_i[nu] + App_i[0] ) <= 0 forall nu > 0
            constraint_name.str(std::string()); // clear the name
            constraint_name << "disappearance/appearance coupling: ";
            constraint_name << " V_i[nu] = 1 => A_i[nu] = 1 v A_i[0] = 1 added for nodes "
                            << constraint.appearance_node << ", " << constraint.disappearance_node;
            constraint_name << " for state: " << state;
            optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), -1,
                    0, constraint_name.str().c_str());
            LOG(logDEBUG3) << constraint_name.str();
        }

        if (!with_misdetections_) {
            cplex_idxs.clear();
            coeffs.clear();

            // assume both nodes are given
            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.disappearance_node, 0));
            coeffs.push_back(1);

            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.appearance_node, 0));
            coeffs.push_back(1);

            // V_i[0] = 0 => 1 <= A_i[0]
            // A_i[0] = 0 => 1 <= V_i[0]
            // V_i <= m, A_i <= m
            // 0 <= Dis_i[0] + App_i[0] <= 0
            constraint_name.str(std::string()); // clear the name
            constraint_name << "disappearance/appearance coupling: ";
            constraint_name << " A_i[0] + V_i[0] = 0 added for nodes " << constraint.appearance_node << ", " << constraint.disappearance_node;
            optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 0,
                    constraint_name.str().c_str());
            LOG(logDEBUG3) << constraint_name.str();
        }

        if (!with_disappearance_) {
            cplex_idxs.clear();
            coeffs.clear();
            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.disappearance_node, 0));
            coeffs.push_back(1);
            // V_i[0] = 0
            // 1 <= V_i <= m
            constraint_name.str(std::string()); // clear the name
            constraint_name << "disappearance/appearance coupling: ";
            constraint_name << " V_i[0] = 0 added for n = "
                    << constraint.disappearance_node;
            optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0,
                    0, constraint_name.str().c_str());
            LOG(logDEBUG3) << constraint_name.str();
        }

        if (!with_appearance_) {
            cplex_idxs.clear();
            coeffs.clear();
            cplex_idxs.push_back(optimizer.lpNodeVi(constraint.appearance_node, 0));
            coeffs.push_back(1);
            // A_i[0] = 0
            // 1 <= A_i <= m
            constraint_name.str(std::string()); // clear the name
            constraint_name << "disappearance/appearance coupling: ";
            constraint_name << " A_i[0] = 0 added for n = "
                    << constraint.appearance_node;
            optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0,
                    0, constraint_name.str().c_str());
            LOG(logDEBUG3) << constraint_name.str();
        }
    }
}

//------------------------------------------------------------------------
template<class CONSTRAINT_TYPE>
void ConstraintPool::constraint_indices(std::vector<ConstraintPool::IndexType>&, const CONSTRAINT_TYPE&)
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

    // if division node id < 0, don't add it and use OutgoingNoDivConstraintFunction
    if(constraint.division_node >= 0)
        indices.push_back(constraint.division_node);

    indices.insert(indices.begin(), constraint.transition_nodes.begin(), constraint.transition_nodes.end());
}

template<>
void ConstraintPool::constraint_indices(std::vector<ConstraintPool::IndexType>& indices, const DetectionConstraint& constraint)
{
    indices.push_back(constraint.disappearance_node);
    indices.push_back(constraint.appearance_node);
}

//------------------------------------------------------------------------
template<class FUNCTION_TYPE>
void ConstraintPool::configure_function(FUNCTION_TYPE* func)
{
    throw std::logic_error("only template specializations of this method should be called");
}

template<>
void ConstraintPool::configure_function(IncomingConstraintFunction<ValueType, IndexType, LabelType>*)
{
    // no flags needed
}

template<>
void ConstraintPool::configure_function(OutgoingNoDivConstraintFunction<ValueType, IndexType, LabelType>*)
{
    // no flags needed
}

template<>
void ConstraintPool::configure_function(OutgoingConstraintFunction<ValueType, IndexType, LabelType>* func)
{
    func->set_with_divisions(with_divisions_);
}

template<>
void ConstraintPool::configure_function(DetectionConstraintFunction<ValueType, IndexType, LabelType>* func)
{
    func->set_with_appearance(with_appearance_);
    func->set_with_disappearance(with_disappearance_);
    func->set_with_misdetections(with_misdetections_);
}

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

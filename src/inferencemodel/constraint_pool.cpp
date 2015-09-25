#include "pgmlink/inferencemodel/constraint_pool.hxx"

namespace pgmlink
{
namespace pgm
{

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
    {
        outgoing_constraints_.push_back(constraint);
    }
    else
    {
        outgoing_no_div_constraints_.push_back(constraint);
    }
}

template<>
void ConstraintPool::add_constraint(const ConstraintPool::DetectionConstraint& constraint)
{
    detection_constraints_.push_back(constraint);
}

template<>
void ConstraintPool::add_constraint(const ConstraintPool::FixNodeValueConstraint& constraint)
{
    fix_node_value_constraints_.push_back(constraint);
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
     )
{
    LOG(logINFO) << "[ConstraintPool]: Adding " << constraints.size() << " hard constraints for Incoming";
    for(auto it = constraints.begin(); it != constraints.end(); ++it)
    {
        const ConstraintPool::IncomingConstraint& constraint = *it;

        // nothing to do if no incoming
        if(constraint.transition_nodes.size() == 0)
        {
            continue;
        }

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
            coeffs.push_back(-(int)state);
        }

        optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 0, 0, constraint_name.str().c_str());
        LOG(logDEBUG3) << constraint_name.str();
    }
}

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
     )
{
    LOG(logINFO) << "[ConstraintPool]: Adding " << constraints.size() << " hard constraints for Outgoing";
    for(auto it = constraints.begin(); it != constraints.end(); ++it)
    {
        const ConstraintPool::OutgoingConstraint& constraint = *it;

        // nothing to do if no outgoing
        if(constraint.transition_nodes.size() == 0)
        {
            continue;
        }

        std::vector<size_t> cplex_idxs, cplex_idxs2;
        std::vector<int> coeffs, coeffs2;

        std::stringstream constraint_name;
        constraint_name << "outgoing transitions: sum( transition-nodes ";

        // couple detection and transitions: Y_ij <= App_i
        for (auto outgoing_it = constraint.transition_nodes.begin(); outgoing_it != constraint.transition_nodes.end(); ++outgoing_it)
        {
            for(size_t a_state = 0; a_state < model.numberOfLabels(constraint.appearance_node); a_state++)
            {
                for(size_t t_state = a_state + 1; t_state < model.numberOfLabels(*outgoing_it) - 1; t_state++)
                {
                    cplex_idxs.clear();
                    coeffs.clear();
                    coeffs.push_back(1);
                    cplex_idxs.push_back(optimizer.lpNodeVi(constraint.appearance_node, a_state));
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
            constraint_name.str(std::string()); // clear the name
            constraint_name << "couple transitions: sum( transition-nodes ";

            for (auto outgoing_it = constraint.transition_nodes.begin(); outgoing_it != constraint.transition_nodes.end(); ++outgoing_it)
            {
                constraint_name << *outgoing_it << " ";
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
                coeffs.push_back(-(int)state);
                cplex_idxs.push_back(optimizer.lpNodeVi(constraint.appearance_node, state));
            }

            // 0 <= sum_nu [ sum_j( nu * Y_ij[nu] ) ] - [ sum_nu nu * X_i[nu] + D_i[1] + sum_nu nu * App_i[nu] ]<= 0
            constraint_name << ") = D_i + App_i added for nodes: App_i=" << constraint.appearance_node
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
                                    -int(model.numberOfLabels(constraint.appearance_node) - 1), 0, constraint_name.str().c_str());
            LOG(logDEBUG3) << constraint_name.str();
        }
    }
}

template<>
void ConstraintPool::add_constraint_type_to_problem<ConstraintPoolOpengmModel,
     ConstraintPoolCplexOptimizer,
     OutgoingNoDivConstraintFunction<ConstraintPool::ValueType, ConstraintPool::IndexType, ConstraintPool::LabelType>,
     ConstraintPool::OutgoingConstraint>
     (
         ConstraintPoolOpengmModel& model,
         ConstraintPoolCplexOptimizer& optimizer,
         const std::vector<ConstraintPool::OutgoingConstraint>& constraints
     )
{
    // for the CPLEX specialization we do the same for with and without divisions
    add_constraint_type_to_problem<ConstraintPoolOpengmModel,
                                   ConstraintPoolCplexOptimizer,
                                   OutgoingConstraintFunction<ValueType, IndexType, LabelType>,
                                   OutgoingConstraint>
                                   (model, optimizer, constraints);
}

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
     )
{
    LOG(logINFO) << "[ConstraintPool]: Adding " << constraints.size() << " hard constraints for Detection";
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

        if (!with_misdetections_)
        {
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

        if (!with_disappearance_)
        {
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

        if (!with_appearance_)
        {
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
// specialization for FixNodeValueConstraintFunction
template<>
void ConstraintPool::add_constraint_type_to_problem<ConstraintPoolOpengmModel,
     ConstraintPoolCplexOptimizer,
     FixNodeValueConstraintFunction<ConstraintPool::ValueType, ConstraintPool::IndexType, ConstraintPool::LabelType>,
     ConstraintPool::FixNodeValueConstraint>
     (
         ConstraintPoolOpengmModel& model,
         ConstraintPoolCplexOptimizer& optimizer,
         const std::vector<ConstraintPool::FixNodeValueConstraint>& constraints
     )
{
    LOG(logINFO) << "[ConstraintPool]: Adding " << constraints.size() << " hard constraints for FixNodeValue";
    for(auto it = constraints.begin(); it != constraints.end(); ++it)
    {
        const ConstraintPool::FixNodeValueConstraint& constraint = *it;

        std::vector<size_t> cplex_idxs;
        std::vector<int> coeffs;
        std::stringstream constraint_name;
        constraint_name << "fix node value of " << constraint.node << " to " << constraint.value;

        cplex_idxs.push_back(optimizer.lpNodeVi(constraint.node, constraint.value));
        coeffs.push_back(1);

        optimizer.addConstraint(cplex_idxs.begin(), cplex_idxs.end(), coeffs.begin(), 1, 1, constraint_name.str().c_str());
        LOG(logDEBUG3) << constraint_name.str();
    }
}

//------------------------------------------------------------------------
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
    {
        indices.push_back(constraint.division_node);
    }

    indices.insert(indices.begin(), constraint.transition_nodes.begin(), constraint.transition_nodes.end());
}

template<>
void ConstraintPool::constraint_indices(std::vector<ConstraintPool::IndexType>& indices, const DetectionConstraint& constraint)
{
    indices.push_back(constraint.disappearance_node);
    indices.push_back(constraint.appearance_node);
}

template<>
void ConstraintPool::constraint_indices(std::vector<ConstraintPool::IndexType>& indices, const FixNodeValueConstraint& constraint)
{
    indices.push_back(constraint.node);
}

//------------------------------------------------------------------------
template<>
void ConstraintPool::configure_function(IncomingConstraintFunction<ValueType, IndexType, LabelType>*, IncomingConstraint)
{
    // no flags needed
}

template<>
void ConstraintPool::configure_function(OutgoingNoDivConstraintFunction<ValueType, IndexType, LabelType>*, OutgoingConstraint)
{
    // no flags needed
}

template<>
void ConstraintPool::configure_function(OutgoingConstraintFunction<ValueType, IndexType, LabelType>* func, IncomingConstraint)
{
    func->set_with_divisions(with_divisions_);
}

template<>
void ConstraintPool::configure_function(DetectionConstraintFunction<ValueType, IndexType, LabelType>* func, DetectionConstraint)
{
    func->set_with_appearance(with_appearance_);
    func->set_with_disappearance(with_disappearance_);
    func->set_with_misdetections(with_misdetections_);
}

template<>
void ConstraintPool::configure_function(FixNodeValueConstraintFunction<ValueType, IndexType, LabelType>* func, ConstraintPool::FixNodeValueConstraint constraint)
{
    func->set_desired_value(constraint.value);
}

} // namespace pgm
} // namespace pgmlink

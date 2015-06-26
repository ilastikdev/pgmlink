//#include "pgmlink/inferencemodel/constrackinginferencemodel.h"
#include <opengm/functions/learnable/lsum_of_experts.hxx>
#include "pgmlink/inferencemodel/structuredlearningtrackinginferencemodel.h"
#include <boost/python.hpp>

//#include "opengm/opengm.hxx"

namespace pgmlink
{

//StructuredLearningTrackingInferenceModel::StructuredLearningTrackingInferenceModel(const Parameter& param,
//        double ep_gap,
//        double cplex_timeout):
//    ConsTrackingInferenceModel(param, ep_gap, cplex_timeout)//,
//    //number_of_transition_nodes_(0),
//    //number_of_division_nodes_(0),
//    //number_of_appearance_nodes_(0),
//    //number_of_disappearance_nodes_(0),
//    //ground_truth_filename_("")
//{
//    //cplex_param_.verbose_ = true;
//    //cplex_param_.integerConstraint_ = true;
//    //cplex_param_.epGap_ = ep_gap;
//    //cplex_param_.timeLimit_ = cplex_timeout;
//}

void StructuredLearningTrackingInferenceModel::setWeight(size_t index, double val){
    weights_[index] = val;
    std::cout << " ===================================================" << weights_[index] << std::endl;
}
/*
void StructuredLearningTrackingInferenceModel::build_from_graph(const HypothesesGraph& hypotheses)
{
    std::cout << "______________________________IN____________________________________build_from_graph!!!" << std::endl;
    LOG(logDEBUG) << "ConsTrackingInferenceModel::formulate: entered";

    LOG(logDEBUG) << "ConsTrackingInferenceModel::formulate: add_transition_nodes";
    add_transition_nodes(hypotheses);
    LOG(logDEBUG) << "ConsTrackingInferenceModel::formulate: add_appearance_nodes";
    add_appearance_nodes(hypotheses);
    LOG(logDEBUG) << "ConsTrackingInferenceModel::formulate: add_disappearance_nodes";
    add_disappearance_nodes(hypotheses);

    LOG(logDEBUG) << "ConsTrackingInferenceModel::formulate: add_division_nodes";
    if (param_.with_divisions)
    {
        add_division_nodes(hypotheses);
    }

    LOG(logINFO) << "number_of_transition_nodes_ = " << number_of_transition_nodes_;
    LOG(logINFO) << "number_of_appearance_nodes_ = " << number_of_appearance_nodes_;
    LOG(logINFO) << "number_of_disappearance_nodes_ = " << number_of_disappearance_nodes_;
    LOG(logINFO) << "number_of_division_nodes_ = " << number_of_division_nodes_;

    add_finite_factors(hypotheses);
    add_constraints_to_pool(hypotheses);
    std::cout << "________________________________OUT__________________________________build_from_graph!!!" << std::endl;
}
*/
size_t StructuredLearningTrackingInferenceModel::add_detection_factors(const HypothesesGraph& g, size_t factorIndex)
{
    std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors" << std::endl;
    ////
    //// add detection factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors" << std::endl;
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());
    std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors" << std::endl;

    LOG(logDEBUG) << "ConsTrackingInferenceModel::add_finite_factors: add detection factors";
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors node" << std::endl;
        size_t num_vars = 0;
        std::vector<size_t> vi;
        std::vector<double> cost;

        int node_begin_time = -1;
        int node_end_time = -1;
        std::cout << ".......0" << std::endl;
        if (param_.with_tracklets)
        {
            std::cout << ".......00" << std::endl;
            node_begin_time = tracklet_map_[n].front().Timestep;
            std::cout << ".......001" << std::endl;
            node_end_time = tracklet_map_[n].back().Timestep;
            std::cout << ".......002" << std::endl;
        }
        else
        {
            std::cout << ".......000" << std::endl;
            node_begin_time = traxel_map_[n].Timestep;
            node_end_time = traxel_map_[n].Timestep;
        }


        std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors node 1" << std::endl;
        double energy, e, f, w;
        if (app_node_map_.count(n) > 0)
        {
            std::cout << ".......1" << std::endl;
            vi.push_back(app_node_map_[n]);
            std::cout << ".......2" << std::endl;
            // "<" holds if there are only tracklets in the first frame
            if (node_begin_time <= g.earliest_timestep())
            {
                std::cout << ".......3" << std::endl;
                // pay no appearance costs in the first timestep
                cost.push_back(0.);
                std::cout << ".......4" << std::endl;
            }
            else
            {
                if (param_.with_tracklets)
                {
                    std::cout << ".......5" << std::endl;
                    energy = param_.appearance_cost(tracklet_map_[n].front());
                    std::cout << ".......6" << std::endl;
                }
                else
                {
                    std::cout << ".......7" << std::endl;
                    energy = param_.appearance_cost(traxel_map_[n]);
                    std::cout << ".......8" << std::endl;
                }

                std::cout << ".......8.0" << std::endl;
                energy += generateRandomOffset(Appearance);
                std::cout << ".......8.1" << std::endl;
                cost.push_back(energy);
                std::cout << ".......8.2" << std::endl;
            }
            std::cout << ".......9" << std::endl;
            ++num_vars;
        }

        std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors node 2" << std::endl;
        if (dis_node_map_.count(n) > 0)
        {
            vi.push_back(dis_node_map_[n]);
            // "<" holds if there are only tracklets in the last frame
            if (node_end_time < g.latest_timestep())
            {
                if (param_.with_tracklets)
                {
                    energy = param_.disappearance_cost(tracklet_map_[n].back());
                }
                else
                {
                    energy = param_.disappearance_cost(traxel_map_[n]);
                }

                energy += generateRandomOffset(Disappearance);
                cost.push_back(energy);
            }
            else
            {
                cost.push_back(0);
            }
            ++num_vars;
        }
        std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors node 3" << std::endl;
        // convert vector to array
        std::vector<size_t> coords(num_vars, 0); // number of variables
        // ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
        //pgm::OpengmExplicitFactor<double> table(vi.begin(), vi.end(),
        // param_.forbidden_cost, (param_.max_number_objects + 1));
        std::vector<size_t> shape(num_vars, (param_.max_number_objects + 1));
        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);
        //marray::Marray<double> energiesP(shape.begin(), shape.end(), param_.forbidden_cost);
        //marray::Marray<double> energiesA(shape.begin(), shape.end(), param_.forbidden_cost);
        //marray::Marray<double> energiesD(shape.begin(), shape.end(), param_.forbidden_cost);
        //marray::Marray<double> energiesT(shape.begin(), shape.end(), param_.forbidden_cost);
        //marray::Marray<double> energiesDiv(shape.begin(), shape.end(), param_.forbidden_cost);

        std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors node 4" << std::endl;
        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
            std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors node state" << std::endl;
            if (param_.with_tracklets)
            {
                energy = 0;
                // add all detection factors of the internal nodes
                for (std::vector<Traxel>::const_iterator trax_it = tracklet_map_[n].begin();
                        trax_it != tracklet_map_[n].end(); ++trax_it)
                {
                    e = param_.detection(*trax_it, state);
                    energy += e;

                    energy += generateRandomOffset(Detection, e, *trax_it, 0, state);
                }
                // add all transition factors of the internal arcs
                Traxel tr_prev;
                bool first = true;
                for (std::vector<Traxel>::const_iterator trax_it = tracklet_map_[n].begin();
                        trax_it != tracklet_map_[n].end(); ++trax_it)
                {
                    LOG(logDEBUG4) << "internal arcs traxel " << *trax_it;
                    Traxel tr = *trax_it;
                    if (!first)
                    {
                        e = param_.transition( get_transition_probability(tr_prev, tr, state) );
                        energy += e;
                        //energy += generateRandomOffset(Transition, e, tr_prev, tr);
                    }
                    else
                    {
                        first = false;
                    }
                    tr_prev = tr;
                }

            }
            else
            {
                e = param_.detection(traxel_map_[n], state);
                std::cout << "detection=e   >" << e << std::endl;
                f = param_.detectionNoWeight(traxel_map_[n], state);
                std::cout << "detection=e/w >" << f << std::endl;

                w = conservationParam_.detection_weight;
                std::cout << "detection=w   >" << w << std::endl;

                energy = e;//f

                //energy += generateRandomOffset(Detection, e, traxel_map_[n], 0, state);
            }
            std::cout << "===>" << energy << std::endl;
            LOG(logDEBUG2) << "ConsTrackingInferenceModel::add_finite_factors: detection[" << state
                           << "] = " << energy;
            for (size_t var_idx = 0; var_idx < num_vars; ++var_idx)
            {
                coords[var_idx] = state;
                // if only one of the variables is > 0, then it is an appearance in this time frame
                // or a disappearance in the next timeframe. Hence, add the cost of appearance/disappearance
                // to the detection cost


                energies(coords.begin()) = conservationParam_.detection_weight * energy + state * cost[var_idx]; // state == m
                //std::cout << "___> " << conservationParam_.detection_weight << "   " << energy << "   " << state << "    " << cost[var_idx] << "   " << var_idx <<std::endl;
                //energiesP(coords.begin()) = energy;




                //if (var_idx == 0) // A
                //    energiesA(coords.begin()) = state;// * cost[var_idx];
                //else // D
                //    energiesD(coords.begin()) = state;// * cost[var_idx];






                coords[var_idx] = 0;
                LOG(logDEBUG4) << "ConsTrackingInferenceModel::add_finite_factors: var_idx "
                               << var_idx << " = " << energy;
            }
            // also this energy if both variables have the same state
            if (num_vars == 2)
            {
                coords[0] = state;
                coords[1] = state;
                // only pay detection energy if both variables are on
                energies(coords.begin()) = conservationParam_.detection_weight * energy;

                //energiesP(coords.begin()) = energy;
                coords[0] = 0;
                coords[1] = 0;

                LOG(logDEBUG4) << "ConsTrackingInferenceModel::add_finite_factors: var_idxs 0 and var_idx 1 = "
                               << energy;
            }
        } // end for state

        std::cout << "energies ==========================================================" << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energies(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << "energiesP ................." << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energiesP(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << "energiesA ................." << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energiesA(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << "energiesD ................." << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energiesD(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }

//        marray::Marray<double> errorMA(shape.begin(), shape.end(), param_.forbidden_cost);
//        errorMA = energies - conservationParam_.detection_weight * energiesP - cost[0]*energiesA - cost[1]*energiesD;



        // a decomposition of energies:
        // w_detection = conservationParam_.detection_weight
        // w_appearance = cost[0]
        // w_disappearance = cost[1]

        // set up a LEARNABLE OpenGM function to learn weights=(w_detection,w_appearance,w_disappearance) for functions (energiesP,energiesA,energiesD)

//        std::cout << "errorMA ................." << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << errorMA(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }


        //for (size_t state = 0; state <= param_.max_number_objects; ++state){
        //    for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
        //        coords[0] = state;
        //        coords[1] = state2;
        //        if (errorMA(coords.begin()) > 0.001 )
        //            std::cout << " We have a problem!!!";
        //    }
        //    //std::cout << std::endl;
        //}

        LOG(logDEBUG3) << "ConsTrackingInferenceModel::add_finite_factors: adding table to pgm";
        //functor add detection table
        factorIndex = add_div_m_best_perturbation(energies, Detection, factorIndex);


//        std::vector<size_t> weightIDs;
//        weightIDs.push_back((size_t)0);
//        weightIDs.push_back((size_t)1);
//        weightIDs.push_back((size_t)2);
//        //for(int i=0;i<3;i++)
//        //    std::cout << weights_[weightIDs[i]] << std::endl;

//        std::vector<marray::Marray<double>> features;
//        features.push_back(energiesP);
//        features.push_back(energiesA);
//        features.push_back(energiesD);
//        //features.push_back(energiesT);
//        //features.push_back(energiesDiv);
//        opengm::functions::learnable::LSumOfExperts<double,size_t,size_t> funEnergies (shape,weights_,weightIDs,features);






        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(energies);
        //typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);

        sort(vi.begin(), vi.end());
        model_.addFactor(funcId, vi.begin(), vi.end());
//        if (not perturb)
        detection_f_node_map_[n] = model_.numberOfFactors() - 1;
    } // end for node n

    return factorIndex;
}

size_t StructuredLearningTrackingInferenceModel::add_transition_factors(const HypothesesGraph& g, size_t factorIndex)
{
    std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_transition_factors" << std::endl;
    ////
    //// add transition factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());

    LOG(logDEBUG) << "ConsTrackingInferenceModel::add_finite_factors: add transition factors";

    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        size_t vi[] = { arc_map_[a] };
        std::vector<size_t> coords(1, 0); // number of variables

        std::vector<size_t> shape(1, (param_.max_number_objects + 1));
        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
            Traxel tr1, tr2;
            if (param_.with_tracklets)
            {
                tr1 = tracklet_map_[g.source(a)].back();
                tr2 = tracklet_map_[g.target(a)].front();
            }
            else
            {
                tr1 = traxel_map_[g.source(a)];
                tr2 = traxel_map_[g.target(a)];
            }

            double energy = param_.transition(get_transition_probability(tr1, tr2, state));
            energy += generateRandomOffset(Transition, energy, tr1, tr2);

            LOG(logDEBUG2) << "ConsTrackingInferenceModel::add_finite_factors: transition[" << state
                           << "] = " << energy;
            coords[0] = state;
            energies(coords.begin()) = energy;
            coords[0] = 0;
        }
        factorIndex = add_div_m_best_perturbation(energies, Transition, factorIndex);


//        std::vector<size_t> weightIDs;
//        weightIDs.push_back((size_t)3);
//        //for(int i=0;i<1;i++)
//        //    std::cout << weights_[weightIDs[i]] << std::endl;
//        std::vector<marray::Marray<double>> features;
//        //features.push_back(energiesP);
//        //features.push_back(energiesA);
//        //features.push_back(energiesD);
//        features.push_back(energies);
//        //features.push_back(energiesDiv);
//        opengm::functions::learnable::LSumOfExperts<double,size_t,size_t> funEnergies (shape,weights_,weightIDs,features);





//        //typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);

        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(energies);
        model_.addFactor(funcId, vi, vi + 1);
    }

    return factorIndex;
}

size_t StructuredLearningTrackingInferenceModel::add_division_factors(const HypothesesGraph& g, size_t factorIndex)
{
    if(!param_.with_divisions)
    {
        return factorIndex;
    }

    std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_division_factors" << std::endl;
    ////
    //// add division factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());

    LOG(logDEBUG) << "ConsTrackingInferenceModel::add_finite_factors: add division factors";
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        if (div_node_map_.count(n) == 0)
        {
            continue;
        }
        size_t vi[] = { div_node_map_[n] };
        std::vector<size_t> coords(1, 0); // number of variables
        std::vector<size_t> shape(1, 2);
        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= 1; ++state)
        {
            double energy;
            Traxel tr;
            if (param_.with_tracklets)
            {
                tr = tracklet_map_[n].back();
            }
            else
            {
                tr = traxel_map_[n];
            }
            energy = param_.division(tr, state);
            energy += generateRandomOffset(Division, energy,  tr);

            LOG(logDEBUG2) << "ConsTrackingInferenceModel::add_finite_factors: division[" << state
                           << "] = " << energy;
            coords[0] = state;
            //table.set_value(coords, energy);
            energies(coords.begin()) = energy;
            coords[0] = 0;
        }
        //table.add_to(model);
        factorIndex = add_div_m_best_perturbation(energies, Division, factorIndex);


//        std::vector<size_t> weightIDs;
//        weightIDs.push_back((size_t)4);
//        //for(int i=0;i<1;i++)
//        //    std::cout << weights_[weightIDs[i]] << std::endl;
//        std::vector<marray::Marray<double>> features;
//        //features.push_back(energiesP);
//        //features.push_back(energiesA);
//        //features.push_back(energiesD);
//        //features.push_back(energiesT);
//        features.push_back(energies);
//        opengm::functions::learnable::LSumOfExperts<double,size_t,size_t> funEnergies (shape,weights_,weightIDs,features);








//        //typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);
        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(energies);
        model_.addFactor(funcId, vi, vi + 1);
    }

    return factorIndex;
}

void StructuredLearningTrackingInferenceModel::add_finite_factors(const HypothesesGraph& g)
{
    std::cout << "-------------------------------StructuredLearningTrackingInferenceModel::add_finite_factors" << std::endl;
    // refactor this:

    // we could use this method to calculate the label-specific offset, also.
    // in order to do so, we would have to write a functor that assigns
    // either the energy or the offset to a table
    // this table is than treated either as an explicit factor table
    // or as an offset marray.
    //
    // Note that this makes it necessary to ensure that factors are added nowhere else in the code
    //
    // Also, this implies that there is some functor choosing either energy or offset

    LOG(logDEBUG) << "ConsTrackingInferenceModel::add_finite_factors: entered";
    size_t factorIndex = 0;
    factorIndex = add_detection_factors(g, factorIndex);
    factorIndex = add_transition_factors(g, factorIndex);
    factorIndex = add_division_factors(g, factorIndex);
    LOG(logDEBUG) << "ConsTrackingInferenceModel::add_finite_factors: finished";
}

//set up optimizer from constraints by reading from formulated gm
void StructuredLearningTrackingInferenceModel::add_constraints_to_pool(const HypothesesGraph& g)
{
    LOG(logDEBUG) << "ConsTrackingInferenceModel::add_constraints: entered";

    constraint_pool_ = pgm::ConstraintPool(param_.forbidden_cost,
                                           param_.with_divisions,
                                           param_.with_appearance,
                                           param_.with_disappearance,
                                           param_.with_misdetections_allowed);

    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = g.get(node_traxel());

    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        std::cout << "add_constraints_to_pool - Time: " << timestep_map[n] << " Node: " << traxel_map[n].Id << std::endl;
        ////
        //// outgoing transitions
        ////
        {
            std::vector<size_t> transition_nodes;
            for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a)
            {
                transition_nodes.push_back(arc_map_[a]);
            }

            int division_node = -1;
            if(div_node_map_.count(n) > 0)
            {
                division_node = div_node_map_[n];
            }
            size_t appearance_node = app_node_map_[n];

            constraint_pool_.add_constraint(pgm::ConstraintPool::OutgoingConstraint(appearance_node,
                                            division_node,
                                            transition_nodes));
        }

        ////
        //// incoming transitions
        ////
        {
            std::vector<size_t> transition_nodes;
            for (HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a)
            {
                transition_nodes.push_back(arc_map_[a]);
            }
            size_t disappearance_node = dis_node_map_[n];

            constraint_pool_.add_constraint(pgm::ConstraintPool::IncomingConstraint(transition_nodes,
                                            disappearance_node));
        }

        ////
        //// disappearance/appearance coupling
        ////
        if (app_node_map_.count(n) > 0 && dis_node_map_.count(n) > 0)
        {
            constraint_pool_.add_constraint(pgm::ConstraintPool::DetectionConstraint((size_t)dis_node_map_[n],
                                            (size_t)app_node_map_[n]));
        }
    }

    constraint_pool_.force_softconstraint(!param_.with_constraints);
}

} // namespace pgmlink

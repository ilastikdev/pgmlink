#ifndef OPENGM_UNSIGNED_INTEGER_POW_HXX_
#define OPENGM_UNSIGNED_INTEGER_POW_HXX_
#endif

#include <opengm/functions/learnable/lweightedsum_of_functions.hxx>
#include "pgmlink/inferencemodel/structuredlearningtrackinginferencemodel.h"
#include <boost/python.hpp>

namespace pgmlink
{

void StructuredLearningTrackingInferenceModel::fixFirstDisappearanceNodesToLabels(
        const HypothesesGraph &g,
        const HypothesesGraph& tracklet_graph,
        std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> >& traxel2tracklet_map
        )
{
    assert(g.has_property(appearance_label()));
    property_map<appearance_label, HypothesesGraph::base_graph>::type &appearance_labels = g.get(appearance_label());

    if(!param_.with_tracklets)
    {
        property_map<node_timestep, HypothesesGraph::base_graph>::type &timestep_map = g.get(node_timestep());
        int earliest_timestep = *(timestep_map.beginValue());

        for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
        {
            if(timestep_map[n] == earliest_timestep)
            {
                linear_constraint_pool_.add_constraint(pgm::ConstraintPool::FixNodeValueLinearConstraint(app_node_map_[n], appearance_labels[n]));
                linear_constraint_pool_.add_constraint(pgm::ConstraintPool::FixNodeValueLinearConstraint(dis_node_map_[n], appearance_labels[n]));
            }
        }
    }
    else
    {
        // in the tracklet graph, the respective label is overwritten by later traxels in the tracklet,
        // get the first original node and use its label
        property_map<node_timestep, HypothesesGraph::base_graph>::type &timestep_map = tracklet_graph.get(node_timestep());
        int earliest_timestep = *(timestep_map.beginValue());

        for (HypothesesGraph::NodeIt n(tracklet_graph); n != lemon::INVALID; ++n)
        {
            if(timestep_map[n] == earliest_timestep)
            {
                HypothesesGraph::Node orig_n = traxel2tracklet_map[n][0];
                linear_constraint_pool_.add_constraint(pgm::ConstraintPool::FixNodeValueLinearConstraint(app_node_map_[n], appearance_labels[orig_n]));
                linear_constraint_pool_.add_constraint(pgm::ConstraintPool::FixNodeValueLinearConstraint(dis_node_map_[n], appearance_labels[orig_n]));
            }
        }
    }
}
size_t StructuredLearningTrackingInferenceModel::add_detection_factors(const HypothesesGraph& g, size_t factorIndex)
{
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());

    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_finite_factors: add detection factors";

    double minEnergyP = 0.0;
    double maxEnergyP = 0.0;
    double minEnergyA = 0.0;
    double maxEnergyA = 0.0;
    double minEnergyD = 0.0;
    double maxEnergyD = 0.0;

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        size_t num_vars = 0;
        double energy;
        if (app_node_map_.count(n) > 0)
        {
            ++num_vars;
        }

        if (dis_node_map_.count(n) > 0)
        {
            ++num_vars;
        }

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
                energy = param_.detectionNoWeight(traxel_map_[n], state);

            for (size_t var_idx = 0; var_idx < num_vars; ++var_idx)
            {
                if (energy < minEnergyP)
                    minEnergyP = energy;
                if (energy > maxEnergyP)
                    maxEnergyP = energy;

                if (var_idx == 0){
                    if (state < minEnergyA)
                        minEnergyA = state;
                    if (state > maxEnergyA)
                        maxEnergyA = state;
                }
                else{
                    if (state < minEnergyD)
                        minEnergyD = state;
                    if (state > maxEnergyD)
                        maxEnergyD = state;
                }
            }
            // also this energy if both variables have the same state
            if (num_vars == 2)
            {
                if (energy < minEnergyP)
                    minEnergyP = energy;
                if (energy > maxEnergyP)
                    maxEnergyP = energy;
            }
        } // end for state
    } // end for node n

    double counterApp=0;
    double counterDis=0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {

        size_t num_vars = 0;
        std::vector<size_t> vi;
        std::vector<double> cost;

        int node_begin_time = -1;
        int node_end_time = -1;
            node_begin_time = traxel_map_[n].Timestep;
            node_end_time = traxel_map_[n].Timestep;

        double energy;
        if (app_node_map_.count(n) > 0)
        {
            vi.push_back(app_node_map_[n]);
            // "<" holds if there are only tracklets in the first frame
            if (node_begin_time <= getModelStartTime())
            {
                // pay no appearance costs in the first timestep
                cost.push_back(0.);
            }
            else
            {
                energy = param_.appearance_cost_fn(traxel_map_[n]);
                cost.push_back(energy);
            }
            ++num_vars;
        }

        if (dis_node_map_.count(n) > 0)
        {
            vi.push_back(dis_node_map_[n]);
            // "<" holds if there are only tracklets in the last frame
            if (node_end_time < getModelEndTime())
            {
                energy = param_.disappearance_cost_fn(traxel_map_[n]);
                cost.push_back(energy);
            }
            else
            {
                cost.push_back(0);
            }
            ++num_vars;
        }
        std::vector<size_t> coords(num_vars, 0);
        std::vector<size_t> shape(num_vars, (param_.max_number_objects + 1));
        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesP(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesA(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesD(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesPnorm(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesAnorm(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesDnorm(shape.begin(), shape.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
            energy = param_.detectionNoWeight(traxel_map_[n], state);

            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: detection[" << state
                           << "] = " << energy;
            for (size_t var_idx = 0; var_idx < num_vars; ++var_idx)
            {
                coords[var_idx] = state;
                // if only one of the variables is > 0, then it is an appearance in this time frame
                // or a disappearance in the next timeframe. Hence, add the cost of appearance/disappearance
                // to the detection cost

                energies(coords.begin()) = learningWeights_.getWeight((size_t)0) * energy + state * cost[var_idx]; // state == m
                assert(maxEnergyP >= energy);
                assert(energy >= minEnergyP);
                energiesPnorm(coords.begin()) = ((double)energy-minEnergyP)/(maxEnergyP-minEnergyP);
                energiesP(coords.begin()) = energy;

                if (var_idx == 0){ // A
                    assert(maxEnergyA >= state);
                    assert(state >= minEnergyA);
                    if (node_begin_time <= getModelStartTime())
                    {
                        energiesAnorm(coords.begin()) = 0.0;
                        energiesA(coords.begin()) = 0.0;
                    }else{
                        double t = traxel_map_[n].Timestep, x = traxel_map_[n].X(), y = traxel_map_[n].Y(), z = traxel_map_[n].Z();
                        double distance = fov_.spatial_distance_to_border(t,x,y,z,false);

                        if (distance < borderWidth_){
                            double ratio = distance/borderWidth_;
                            energiesAnorm(coords.begin()) = ((double)ratio*state-minEnergyA)/(maxEnergyA-minEnergyA);
                            energiesA(coords.begin()) = ratio*state;
                            if(state==0)
                                ++counterApp;
                        }
                        else{
                            energiesAnorm(coords.begin()) = ((double)state-minEnergyA)/(maxEnergyA-minEnergyA);
                            energiesA(coords.begin()) = state;
                        }
                    }
                }
                else{ // D
                    assert(maxEnergyD >= state);
                    assert(state >= minEnergyD);
                    if (node_end_time < getModelEndTime()){
                        double t = traxel_map_[n].Timestep, x = traxel_map_[n].X(), y = traxel_map_[n].Y(), z = traxel_map_[n].Z();
                        double distance = fov_.spatial_distance_to_border(t,x,y,z,false);

                        if (distance < borderWidth_){
                            double ratio = distance/borderWidth_;
                            energiesDnorm(coords.begin()) = ((double)ratio*state-minEnergyD)/(maxEnergyD-minEnergyD);
                            energiesD(coords.begin()) = ratio*state;
                            if(state==0)
                                ++counterDis;
                        }
                        else{
                            energiesDnorm(coords.begin()) = ((double)state-minEnergyD)/(maxEnergyD-minEnergyD);
                            energiesD(coords.begin()) = state;
                        }
                    }
                    else{
                        energiesDnorm(coords.begin()) = 0.0;
                        energiesD(coords.begin()) = 0.0;
                    }
                }

                coords[var_idx] = 0;
                LOG(logDEBUG4) << "StructuredLearningTrackingInferenceModel::add_finite_factors: var_idx "
                               << var_idx << " = " << energy;
            }
            // also this energy if both variables have the same state
            if (num_vars == 2)
            {
                coords[0] = state;
                coords[1] = state;

                energies(coords.begin()) = learningWeights_.getWeight((size_t)0) * energy;
                assert(maxEnergyP >= energy);
                assert(energy >= minEnergyP);
                energiesPnorm(coords.begin()) = ((double)energy-minEnergyP)/(maxEnergyP-minEnergyP);
                energiesP(coords.begin()) = energy;

                coords[0] = 0;
                coords[1] = 0;

                LOG(logDEBUG4) << "StructuredLearningTrackingInferenceModel::add_finite_factors: var_idxs 0 and var_idx 1 = "
                               << energy;
            }
        } // end for state

        // a decomposition of energies:
        // w_detection = conservationParam_.detection_weight
        // w_appearance = cost[0]
        // w_disappearance = cost[1]

        LOG(logDEBUG3) << "StructuredLearningTrackingInferenceModel::add_finite_factors: adding table to pgm";

        std::vector<size_t> weightIDs;
        weightIDs.push_back((size_t)0);
        weightIDs.push_back((size_t)3);
        weightIDs.push_back((size_t)4);

        std::vector<marray::Marray<double>> featuresNorm;
        featuresNorm.push_back(energiesPnorm);
        featuresNorm.push_back(energiesAnorm);
        featuresNorm.push_back(energiesDnorm);

        std::vector<marray::Marray<double>> features;
        features.push_back(energiesP);
        features.push_back(energiesA);
        features.push_back(energiesD);

        std::vector<size_t> varShape;
        varShape.push_back((size_t)param_.max_number_objects+1);
        varShape.push_back((size_t)param_.max_number_objects+1);

        opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t> funEnergies;
        if(withNormalization_)
            funEnergies = opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t>(varShape,learningWeights_,weightIDs,featuresNorm);
        else
            funEnergies = opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t>(varShape,learningWeights_,weightIDs,features);

        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);

        // sorting only works because appearance nodes have lower variable indices than disappearance nodes
        // and the matrix is constructed such that appearances are along coords[0], ...
        sort(vi.begin(), vi.end());
        model_.addFactor(funcId, vi.begin(), vi.end());
        detection_f_node_map_[n] = model_.numberOfFactors() - 1;
    } // end for node n
    LOG(logDEBUG) << "number of REDUCED border appearance weights    " << counterApp << std::endl;
    LOG(logDEBUG) << "number of REDUCED border disappearance weights " << counterDis << std::endl;
    return factorIndex;
}

size_t StructuredLearningTrackingInferenceModel::add_transition_factors(const HypothesesGraph& g, size_t factorIndex)
{
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());

    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_finite_factors: add transition factors";

    double minEnergy = (double)std::numeric_limits<double>::infinity();
    double maxEnergy = 0.0;

    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        Traxel tr1, tr2;
        tr1 = traxel_map_[g.source(a)];
        tr2 = traxel_map_[g.target(a)];

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
            double energy = param_.transitionNoWeight(get_transition_probability(tr1, tr2, state));
            //double energy = param_.transition(tr1, tr2, state);

            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: transition[" << state
                           << "] = " << energy;
            if (energy < minEnergy)
                minEnergy = energy;
            if (energy > maxEnergy)
                maxEnergy = energy;
        }
    }

    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        size_t vi[] = { arc_map_[a] };
        std::vector<size_t> coords(1, 0);

        std::vector<size_t> shape(1, (param_.max_number_objects + 1));
        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesNorm(shape.begin(), shape.end(), param_.forbidden_cost);

        Traxel tr1, tr2;
        tr1 = traxel_map_[g.source(a)];
        tr2 = traxel_map_[g.target(a)];

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
            double energy = param_.transitionNoWeight(get_transition_probability(tr1, tr2, state));
            //double energy = param_.transition(tr1, tr2, state); // calls python side defined transition funtion WITHOUT weight multiplier

            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: transition[" << state
                           << "] = " << energy;
            coords[0] = state;
            assert(maxEnergy >= energy);
            assert(energy >= minEnergy);
            energiesNorm(coords.begin()) = ((double)energy-minEnergy)/(maxEnergy-minEnergy);
            energies(coords.begin()) = energy;

            coords[0] = 0;
        }
        std::vector<size_t> weightIDs;
        weightIDs.push_back((size_t)2);

        std::vector<marray::Marray<double>> featuresNorm;
        featuresNorm.push_back(energiesNorm);

        std::vector<marray::Marray<double>> features;
        features.push_back(energies);

        std::vector<size_t> varShape;
        varShape.push_back((size_t)1+param_.max_number_objects);

        opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t> funEnergies;
        if(withNormalization_)
            funEnergies = opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t>(varShape,learningWeights_,weightIDs,features);
        else
            funEnergies = opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t>(varShape,learningWeights_,weightIDs,featuresNorm);

        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);
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

    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());

    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_finite_factors: add division factors";

    double minEnergy = (double)std::numeric_limits<double>::infinity();
    double maxEnergy = 0.0;

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        if (div_node_map_.count(n) == 0)
        {
            continue;
        }
        Traxel tr;
        tr = traxel_map_[n];
        for (size_t state = 0; state <= 1; ++state)
        {
            double energy;
            energy = param_.divisionNoWeight(tr, state);

            if (energy < minEnergy)
                minEnergy = energy;
            if (energy > maxEnergy)
                maxEnergy = energy;
        }
    }

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        if (div_node_map_.count(n) == 0)
        {
            continue;
        }
        size_t vi[] = { div_node_map_[n] };
        std::vector<size_t> coords(1, 0);
        std::vector<size_t> shape(1, 2);
        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesNorm(shape.begin(), shape.end(), param_.forbidden_cost);

        Traxel tr;
        tr = traxel_map_[n];

        for (size_t state = 0; state <= 1; ++state)
        {
            double energy;
            energy = param_.divisionNoWeight(tr, state);

            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: division[" << state
                           << "] = " << energy;
            coords[0] = state;
            assert(maxEnergy >= energy);
            assert(energy >= minEnergy);
            energiesNorm(coords.begin()) = ((double)energy-minEnergy)/(maxEnergy-minEnergy);
            energies(coords.begin()) = energy;

            coords[0] = 0;
        }

        std::vector<size_t> weightIDs;
        weightIDs.push_back((size_t)1);

        std::vector<marray::Marray<double>> featuresNorm;
        featuresNorm.push_back(energiesNorm);

        std::vector<marray::Marray<double>> features;
        features.push_back(energies);

        std::vector<size_t> varShape;
        varShape.push_back((size_t)2);

        opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t> funEnergies;
        if(withNormalization_)
            funEnergies = opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t>(varShape,learningWeights_,weightIDs,features);
        else
            funEnergies = opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t>(varShape,learningWeights_,weightIDs,featuresNorm);

        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);
        model_.addFactor(funcId, vi, vi + 1);

    }

    return factorIndex;
}

void StructuredLearningTrackingInferenceModel::add_constraints_to_pool(const HypothesesGraph& g)
{
    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_constraints_to_pool";

    linear_constraint_pool_ = pgm::ConstraintPool(param_.forbidden_cost,
                                           param_.with_divisions,
                                           param_.with_appearance,
                                           param_.with_disappearance,
                                           param_.with_misdetections_allowed);

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
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

            linear_constraint_pool_.add_constraint(pgm::ConstraintPool::OutgoingLinearConstraint(appearance_node,division_node,transition_nodes));
        }

        {
            std::vector<size_t> transition_nodes;
            for (HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a)
            {
                transition_nodes.push_back(arc_map_[a]);
            }
            size_t disappearance_node = dis_node_map_[n];

            linear_constraint_pool_.add_constraint(pgm::ConstraintPool::IncomingLinearConstraint(transition_nodes,disappearance_node));
        }

        if (app_node_map_.count(n) > 0 && dis_node_map_.count(n) > 0)
        {
            linear_constraint_pool_.add_constraint(pgm::ConstraintPool::DetectionLinearConstraint((size_t)dis_node_map_[n],(size_t)app_node_map_[n]));
        }
    }
}

void StructuredLearningTrackingInferenceModel::set_inference_params(size_t numberOfSolutions,
        const std::string &feature_filename,
        const std::string &constraints_filename,
        const std::string &ground_truth_filename)
{
    ground_truth_filename_ = ground_truth_filename;

#ifdef WITH_MODIFIED_OPENGM
    optimizer2_ = boost::shared_ptr<cplex2_optimizer>(new cplex2_optimizer(get_model(),
                 cplex2_param_,
                 numberOfSolutions,
                 feature_filename,
                 constraints_filename,
                 ground_truth_filename));

    cplex_variable_id_map_ = optimizer_->get_cplex_variable_id_map();
    cplex_factor_id_map_ = optimizer_->get_cplex_factor_id_map();
#else
    optimizer2_ = boost::shared_ptr<cplex2_optimizer>(new cplex2_optimizer(get_model(), cplex2_param_));
#endif

    if(param_.with_constraints)
    {
        LOG(logINFO) << "add_constraints";
        add_constraints(*optimizer2_);
    }
    else
    {
        throw std::runtime_error("GraphicalModel::infer(): inference with soft constraints is not implemented yet. "
                                 "The conservation tracking factor graph has been saved to file");
    }
}

ConsTrackingInferenceModel::IlpSolution StructuredLearningTrackingInferenceModel::infer()
{

    opengm::InferenceTermination status = optimizer2_->infer();
    if (status != opengm::NORMAL)
    {
        throw std::runtime_error("GraphicalModel::infer(): optimizer terminated abnormally");
    }

    IlpSolution solution;
    opengm::InferenceTermination statusExtract = optimizer2_->arg(solution);
    if (statusExtract != opengm::NORMAL)
    {
        throw std::runtime_error("GraphicalModel::infer(): solution extraction terminated abnormally");
    }

    return solution;
}

ConsTrackingInferenceModel::IlpSolution StructuredLearningTrackingInferenceModel::extractSolution(size_t k,
        const std::string &ground_truth_filename)
{
    IlpSolution solution;
#ifdef WITH_MODIFIED_OPENGM
    optimizer_->set_export_file_names("", "", ground_truth_filename);
#endif
    opengm::InferenceTermination status = optimizer2_->arg(solution, k);

    if (status != opengm::NORMAL)
    {
        throw std::runtime_error("GraphicalModel::infer(): solution extraction terminated abnormally");
    }

    return solution;
}

} // namespace pgmlink

#include "pgmlink/constrackinginferencemodel.h"
#include <boost/python.hpp>

namespace pgmlink
{

ConsTrackingInferenceModel::ConsTrackingInferenceModel(const Parameter& param,
                                                       double ep_gap,
                                                       double cplex_timeout):
    param_(param),
    number_of_transition_nodes_(0),
    number_of_division_nodes_(0),
    number_of_appearance_nodes_(0),
    number_of_disappearance_nodes_(0),
    transition_predictions_(new TransitionPredictionsMap())
{
    cplex_param_.verbose_ = true;
    cplex_param_.integerConstraint_ = true;
    cplex_param_.epGap_ = ep_gap;
    cplex_param_.timeLimit_ = cplex_timeout;
}

void ConsTrackingInferenceModel::use_transition_prediction_cache(ConsTrackingInferenceModel *other)
{
    if(other == NULL)
    {
        LOG(logWARNING) << "[ConsTrackingInferenceModel] Cannot copy cache from NULL pointer" << std::endl;
    }
    else if(!other->transition_predictions_)
    {
        LOG(logWARNING) << "[ConsTrackingInferenceModel] Not setting empty cache" << std::endl;
    }
    else
    {
        transition_predictions_ = other->transition_predictions_;
    }
}

void ConsTrackingInferenceModel::build_from_graph(const HypothesesGraph& hypotheses) {
    LOG(logDEBUG) << "ConsTrackingInferenceModel::formulate: entered";

    LOG(logDEBUG) << "ConsTrackingInferenceModel::formulate: add_transition_nodes";
    add_transition_nodes(hypotheses);
    LOG(logDEBUG) << "ConsTrackingInferenceModel::formulate: add_appearance_nodes";
    add_appearance_nodes(hypotheses);
    LOG(logDEBUG) << "ConsTrackingInferenceModel::formulate: add_disappearance_nodes";
    add_disappearance_nodes(hypotheses);

    LOG(logDEBUG) << "ConsTrackingInferenceModel::formulate: add_division_nodes";
    if (param_.with_divisions) {
        add_division_nodes(hypotheses);
    }

    LOG(logINFO) << "number_of_transition_nodes_ = " << number_of_transition_nodes_;
    LOG(logINFO) << "number_of_appearance_nodes_ = " << number_of_appearance_nodes_;
    LOG(logINFO) << "number_of_disappearance_nodes_ = " << number_of_disappearance_nodes_;
    LOG(logINFO) << "number_of_division_nodes_ = " << number_of_division_nodes_;

    add_finite_factors(hypotheses);
    add_constraints_to_pool(hypotheses);
}

ConsTrackingInferenceModel::GraphicalModelType& ConsTrackingInferenceModel::get_model()
{
    return model_;
}

ConsTrackingInferenceModel::HypothesesGraphNodeMap &ConsTrackingInferenceModel::get_division_node_map()
{
    return div_node_map_;
}

ConsTrackingInferenceModel::HypothesesGraphNodeMap &ConsTrackingInferenceModel::get_appearance_node_map()
{
    return app_node_map_;
}

ConsTrackingInferenceModel::HypothesesGraphNodeMap &ConsTrackingInferenceModel::get_disappearance_node_map()
{
    return dis_node_map_;
}

ConsTrackingInferenceModel::HypothesesGraphArcMap &ConsTrackingInferenceModel::get_arc_map()
{
    return arc_map_;
}

std::map<HypothesesGraph::Node, size_t>& ConsTrackingInferenceModel::get_detection_factor_node_map()
{
    return detection_f_node_map_;
}

void ConsTrackingInferenceModel::add_appearance_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        model_.addVariable(param_.max_number_objects + 1);
        app_node_map_[n] = model_.numberOfVariables() - 1;

        HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
        nodes_per_timestep_[timestep_map[n]].push_back(model_.numberOfVariables() - 1);

        assert(model_.numberOfLabels(app_node_map_[n]) == param_.max_number_objects + 1);
        ++count;
    }
    number_of_appearance_nodes_ = count;
}

void ConsTrackingInferenceModel::add_disappearance_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        model_.addVariable(param_.max_number_objects + 1);
        dis_node_map_[n] = model_.numberOfVariables() - 1;

        HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
        nodes_per_timestep_[timestep_map[n]].push_back(model_.numberOfVariables() - 1);

        assert(model_.numberOfLabels(dis_node_map_[n]) == param_.max_number_objects + 1);
        ++count;
    }
    number_of_disappearance_nodes_ = count;
}

void ConsTrackingInferenceModel::add_transition_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        model_.addVariable(param_.max_number_objects + 1);
        arc_map_[a] = model_.numberOfVariables() - 1;
        // store these nodes by the timestep of the base-appearance node,
        // as well as in the timestep of the disappearance node they enter
        HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
        HypothesesGraph::OutArcIt out_arc(g, a);
        HypothesesGraph::Node n = g.source(out_arc);
        nodes_per_timestep_[timestep_map[n]].push_back(model_.numberOfVariables() - 1);
        n = g.target(out_arc);
        nodes_per_timestep_[timestep_map[n]].push_back(model_.numberOfVariables() - 1);

        assert(model_.numberOfLabels(arc_map_[a]) == param_.max_number_objects + 1);
        ++count;
    }
    number_of_transition_nodes_ = count;
}

void ConsTrackingInferenceModel::add_division_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        size_t number_of_outarcs = 0;
        for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
            ++number_of_outarcs;
        }
        if (number_of_outarcs > 1) {
            model_.addVariable(2);
            div_node_map_[n] = model_.numberOfVariables() - 1;
            HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
            nodes_per_timestep_[timestep_map[n]].push_back(model_.numberOfVariables() - 1);

            assert(model_.numberOfLabels(div_node_map_[n]) == 2);
            ++count;
        }
    }
    number_of_division_nodes_ = count;
}


void ConsTrackingInferenceModel::printResults(const HypothesesGraph& g)
{
    //very verbose print of solution
    property_map<arc_active_count, HypothesesGraph::base_graph>::type& active_arcs_count =
            g.get(arc_active_count());
    property_map<node_active_count, HypothesesGraph::base_graph>::type& active_nodes_count =
            g.get(node_active_count());
    property_map<division_active_count, HypothesesGraph::base_graph>::type& active_divisions_count =
            g.get(division_active_count());

    int c=0;
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        c=0;
        for( std::vector<bool>::const_iterator i = active_arcs_count[a].begin();
             i != active_arcs_count[a].end();
             ++i){
            LOG(logDEBUG4) << *i;
            if (*i) {
                c++;
            }
        }
        LOG(logDEBUG4) << "total= "<<c;
    }

    for(std::map<HypothesesGraph::Node, size_t>::const_iterator it = app_node_map_.begin();
        it != app_node_map_.end(); ++it) {
        c=0;
        for( std::vector<long unsigned int>::const_iterator i = active_nodes_count[it->first].begin();
             i != active_nodes_count[it->first].end();
             ++i){
            LOG(logINFO) << *i;
            if (*i>0){
                c++;
            }
        }
        LOG(logINFO) << "total= "<<c<<std::endl;
    }
    LOG(logDEBUG4) << "division nodes "<<c<<std::endl;
    for(std::map<HypothesesGraph::Node, size_t>::const_iterator it = app_node_map_.begin();
        it != app_node_map_.end(); ++it) {
        c=0;
        for( std::vector<bool>::const_iterator i = active_divisions_count[it->first].begin();
             i != active_divisions_count[it->first].end();
             ++i){
            LOG(logDEBUG4) << *i<<" ";
            if (*i>0){
                c++;
            }
        }
        LOG(logDEBUG4) << "total= "<<c<<std::endl;
    }
}

double ConsTrackingInferenceModel::get_transition_prob(double distance, size_t state, double alpha) {
    double prob = exp(-distance / alpha);
    if (state == 0) {
        return 1 - prob;
    }
    return prob;
}

double ConsTrackingInferenceModel::generateRandomOffset(EnergyType energyIndex,
                                                        double energy,
                                                        Traxel tr,
                                                        Traxel tr2,
                                                        size_t state)
{
    // unperturbed inference -> no random offset
    return 0.0;
}

bool ConsTrackingInferenceModel::callable(boost::python::object object)
{
  return 1 == PyCallable_Check(object.ptr());
}

double ConsTrackingInferenceModel::get_transition_probability(Traxel& tr1, Traxel& tr2, size_t state) {
    LOG(logDEBUG4) << "get_transition_probability()";

    double prob;

    //read the FeatureMaps from Traxels
    if (param_.transition_classifier.ptr()==boost::python::object().ptr()) {
        double distance = 0;
        if (param_.with_optical_correction) {
            distance = tr1.distance_to_corr(tr2);
        } else {
            distance = tr1.distance_to(tr2);
        }
        prob = get_transition_prob(distance, state, param_.transition_parameter);
        LOG(logDEBUG4) << "get_transition_probability(): using deterministic function: " << tr1
                       << " " << tr2 << " [" << state << "] = " << prob << "; distance = " << distance;
        assert(prob >= 0 && prob <= 1);
        return prob;
    }

    TransitionPredictionsMap::const_iterator it = transition_predictions_->find(std::make_pair(tr1, tr2));
    if ( it == transition_predictions_->end() )
    {
        // predict and store
        double var;
        try {
            assert(tr1.features.find("com") != tr1.features.end());
            assert(callable(param_.transition_classifier.attr("predictWithCoordinates")));
            PyGILState_STATE pygilstate = PyGILState_Ensure();
            boost::python::object prediction =
                    param_.transition_classifier.attr("predictWithCoordinates")(tr1.X(), tr1.Y(), tr1.Z(),
                                                                                tr2.X(), tr2.Y(), tr2.Z());
            boost::python::object probs_python = prediction.attr("__getitem__")(0);
            // we are only interested in the probability of the second class, since it is a binary classifier
            prob = boost::python::extract<double>(probs_python.attr("__getitem__")(1));
            var = boost::python::extract<double>(prediction.attr("__getitem__")(1));
            PyGILState_Release(pygilstate);
        } catch (...) {
            throw std::runtime_error("cannot call the transition classifier from python");
        }
        (*transition_predictions_)[std::make_pair(tr1, tr2)] = std::make_pair(prob, var);
    }

    if (state == 0) {
        prob = 1-prob;
    }
    LOG(logDEBUG4) << "get_transition_probability(): using Gaussian process classifier: p["
                   << state << "] = " << prob ;
    return prob;
}


size_t ConsTrackingInferenceModel::add_detection_factors(const HypothesesGraph& g, size_t factorIndex)
{
    ////
    //// add detection factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
            g.get(node_tracklet());

    LOG(logDEBUG) << "ConsTrackingInferenceModel::add_finite_factors: add detection factors";
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        size_t num_vars = 0;
        std::vector<size_t> vi;
        std::vector<double> cost;

        int node_begin_time = -1;
        int node_end_time = -1;
        if (param_.with_tracklets) {
            node_begin_time = tracklet_map_[n].front().Timestep;
            node_end_time = tracklet_map_[n].back().Timestep;
        } else {
            node_begin_time = traxel_map_[n].Timestep;
            node_end_time = traxel_map_[n].Timestep;
        }


        double energy,e;
        if (app_node_map_.count(n) > 0) {
            vi.push_back(app_node_map_[n]);
            // "<" holds if there are only tracklets in the first frame
            if (node_begin_time <= g.earliest_timestep()) {
                // pay no appearance costs in the first timestep
                cost.push_back(0.);
            } else {
                if (param_.with_tracklets) {
                    energy = param_.appearance_cost(tracklet_map_[n].front());
                } else {
                    energy = param_.appearance_cost(traxel_map_[n]);
                }

                energy += generateRandomOffset(Appearance);
                cost.push_back(energy);
            }
            ++num_vars;
        }

        if (dis_node_map_.count(n) > 0) {
            vi.push_back(dis_node_map_[n]);
            // "<" holds if there are only tracklets in the last frame
            if (node_end_time < g.latest_timestep()) {
                if (param_.with_tracklets) {
                    energy = param_.disappearance_cost(tracklet_map_[n].back());
                } else {
                    energy = param_.disappearance_cost(traxel_map_[n]);
                }

                energy += generateRandomOffset(Disappearance);
                cost.push_back(energy);
            } else {
                cost.push_back(0);
            }
            ++num_vars;
        }
        // convert vector to array
        std::vector<size_t> coords(num_vars, 0); // number of variables
        // ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
        //pgm::OpengmExplicitFactor<double> table(vi.begin(), vi.end(),
        // param_.forbidden_cost, (param_.max_number_objects + 1));
        std::vector<size_t> shape(num_vars,(param_.max_number_objects+1));
        marray::Marray<double> energies(shape.begin(),shape.end(),param_.forbidden_cost);

        for (size_t state = 0; state <= param_.max_number_objects; ++state) {
            if (param_.with_tracklets) {
                energy=0;
                // add all detection factors of the internal nodes
                for (std::vector<Traxel>::const_iterator trax_it = tracklet_map_[n].begin();
                        trax_it != tracklet_map_[n].end(); ++trax_it) {
                    e = param_.detection(*trax_it, state);
                    energy += e;

                    energy += generateRandomOffset(Detection, e, *trax_it, 0, state);
                }
                // add all transition factors of the internal arcs
                Traxel tr_prev;
                bool first = true;
                for (std::vector<Traxel>::const_iterator trax_it = tracklet_map_[n].begin();
                        trax_it != tracklet_map_[n].end(); ++trax_it) {
                    LOG(logDEBUG4) << "internal arcs traxel " << *trax_it;
                    Traxel tr = *trax_it;
                    if (!first) {
                        e = param_.transition( get_transition_probability(tr_prev, tr, state) );
                        energy += e;
                        energy += generateRandomOffset(Transition, e, tr_prev, tr);
                    } else {
                        first = false;
                    }
                    tr_prev = tr;
                }

            } else {
                e = param_.detection(traxel_map_[n], state);
                energy = e;
                energy += generateRandomOffset(Detection, e, traxel_map_[n], 0, state);
            }
            LOG(logDEBUG2) << "ConsTrackingInferenceModel::add_finite_factors: detection[" << state
                 << "] = " << energy;
            for (size_t var_idx = 0; var_idx < num_vars; ++var_idx) {
                coords[var_idx] = state;
                // if only one of the variables is > 0, then it is an appearance in this time frame
                // or a disappearance in the next timeframe. Hence, add the cost of appearance/disappearance
                // to the detection cost
                energies(coords.begin()) = energy + state * cost[var_idx];
                coords[var_idx] = 0;
                LOG(logDEBUG4) << "ConsTrackingInferenceModel::add_finite_factors: var_idx "
                                << var_idx << " = " << energy;
            }
            // also this energy if both variables have the same state
            if (num_vars == 2) {
                coords[0] = state;
                coords[1] = state;
                // only pay detection energy if both variables are on
                energies(coords.begin()) = energy;
                coords[0] = 0;
                coords[1] = 0;

                LOG(logDEBUG4) << "ConsTrackingInferenceModel::add_finite_factors: var_idxs 0 and var_idx 1 = "
                                << energy;
            }
        }

        LOG(logDEBUG3) << "ConsTrackingInferenceModel::add_finite_factors: adding table to pgm";
        //functor add detection table
        factorIndex = add_div_m_best_perturbation(energies, Detection, factorIndex);
        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(energies);

        sort(vi.begin(),vi.end());
        model_.addFactor(funcId,vi.begin(),vi.end());
//        if (not perturb)
        detection_f_node_map_[n] = model_.numberOfFactors() - 1;
    }

    return factorIndex;
}

size_t ConsTrackingInferenceModel::add_transition_factors(const HypothesesGraph& g, size_t factorIndex)
{
    ////
    //// add transition factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
            g.get(node_tracklet());

    LOG(logDEBUG) << "ConsTrackingInferenceModel::add_finite_factors: add transition factors";

    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        size_t vi[] = { arc_map_[a] };
        std::vector<size_t> coords(1, 0); // number of variables

        std::vector<size_t> shape(1,(param_.max_number_objects+1));
        marray::Marray<double> energies(shape.begin(),shape.end(),param_.forbidden_cost);

        for (size_t state = 0; state <= param_.max_number_objects; ++state) {
            Traxel tr1, tr2;
            if (param_.with_tracklets) {
                tr1 = tracklet_map_[g.source(a)].back();
                tr2 = tracklet_map_[g.target(a)].front();
            } else {
                tr1 = traxel_map_[g.source(a)];
                tr2 = traxel_map_[g.target(a)];
            }

            double energy = param_.transition(get_transition_probability(tr1, tr2, state));
            energy += generateRandomOffset(Transition, energy, tr1, tr2);

            LOG(logDEBUG2) << "ConsTrackingInferenceModel::add_finite_factors: transition[" << state
                    << "] = " << energy;
            coords[0] = state;
            energies(coords.begin())=energy;
            coords[0] = 0;
        }
        factorIndex = add_div_m_best_perturbation(energies, Transition, factorIndex);
        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(energies);
        model_.addFactor(funcId,vi,vi+1);
    }

    return factorIndex;
}

size_t ConsTrackingInferenceModel::add_div_m_best_perturbation(marray::Marray<double>& energies,
                                                               EnergyType energy_type,
                                                               size_t factorIndex)
{
    return factorIndex;
}

size_t ConsTrackingInferenceModel::add_division_factors(const HypothesesGraph& g, size_t factorIndex)
{
    if(!param_.with_divisions)
        return factorIndex;

    ////
    //// add division factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
            g.get(node_tracklet());

    LOG(logDEBUG) << "ConsTrackingInferenceModel::add_finite_factors: add division factors";
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        if (div_node_map_.count(n) == 0) {
            continue;
        }
        size_t vi[] = { div_node_map_[n] };
        std::vector<size_t> coords(1, 0); // number of variables
        std::vector<size_t> shape(1,2);
        marray::Marray<double> energies(shape.begin(),shape.end(),param_.forbidden_cost);

        for (size_t state = 0; state <= 1; ++state) {
            double energy;
            Traxel tr;
            if (param_.with_tracklets) {
                tr = tracklet_map_[n].back();
            } else {
                tr = traxel_map_[n];
            }
            energy = param_.division(tr, state);
            energy += generateRandomOffset(Division, energy,  tr);

            LOG(logDEBUG2) << "ConsTrackingInferenceModel::add_finite_factors: division[" << state
                    << "] = " << energy;
            coords[0] = state;
            //table.set_value(coords, energy);
            energies(coords.begin())=energy;
            coords[0] = 0;
        }
        //table.add_to(model);
        factorIndex = add_div_m_best_perturbation(energies, Division, factorIndex);

        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(energies);
        model_.addFactor(funcId,vi,vi+1);
    }

    return factorIndex;
}

void ConsTrackingInferenceModel::add_finite_factors(const HypothesesGraph& g) {
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
void ConsTrackingInferenceModel::add_constraints_to_pool(const HypothesesGraph& g) {
    LOG(logDEBUG) << "ConsTrackingInferenceModel::add_constraints: entered";

    constraint_pool_ = pgm::ConstraintPool(param_.forbidden_cost,
                                           param_.with_divisions,
                                           param_.with_appearance,
                                           param_.with_disappearance,
                                           param_.with_misdetections_allowed);

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        ////
        //// outgoing transitions
        ////
        {
            std::vector<size_t> transition_nodes;
            for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
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
            for (HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a) {
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

ConsTrackingInferenceModel::IlpSolution ConsTrackingInferenceModel::infer(size_t numberOfSolutions,
                                                                          const std::string &feature_filename,
                                                                          const std::string &constraints_filename,
                                                                          const std::string &ground_truth_filename,
                                                                          bool with_inference,
                                                                          bool export_from_labeled_graph)
{
    optimizer_ = boost::shared_ptr<cplex_optimizer>(
            new cplex_optimizer(get_model(),
                                cplex_param_,
                                numberOfSolutions,
                                feature_filename,
                                constraints_filename,
                                ground_truth_filename,
                                export_from_labeled_graph));
    if(!with_inference)
        return IlpSolution();

    if(param_.with_constraints)
    {
        LOG(logINFO) << "add_constraints";
        add_constraints(*optimizer_);
    }
    else
    {
        //opengm::hdf5::save(optimizer_->graphicalModel(), "./conservationTracking.h5", "conservationTracking");
        throw std::runtime_error("GraphicalModel::infer(): inference with soft constraints is not implemented yet. "
                                         "The conservation tracking factor graph has been saved to file");
    }

    opengm::InferenceTermination status = optimizer_->infer();
    if (status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): optimizer terminated abnormally");
    }

    IlpSolution solution;
    opengm::InferenceTermination statusExtract = optimizer_->arg(solution);
    if (statusExtract != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): solution extraction terminated abnormally");
    }
    if(export_from_labeled_graph and not  ground_truth_filename.empty()){
        clpex_variable_id_map_ = optimizer_->get_clpex_variable_id_map();
        clpex_factor_id_map_ = optimizer_->get_clpex_factor_id_map();
    }

    return solution;
}

ConsTrackingInferenceModel::IlpSolution ConsTrackingInferenceModel::extractSolution(size_t k,
                                                        const std::string &ground_truth_filename)
{
    IlpSolution solution;
    optimizer_->set_export_file_names("", "", ground_truth_filename);
    opengm::InferenceTermination status = optimizer_->arg(solution, k);

    if (status != opengm::NORMAL)
    {
        throw std::runtime_error("GraphicalModel::infer(): solution extraction terminated abnormally");
    }

    return solution;
}

void ConsTrackingInferenceModel::write_labeledgraph_to_file(const HypothesesGraph &g,
                                                            const std::string &ground_truth_filename)
{
    using namespace std;

    cout << "write_labeledgraph_to_file" << endl;
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    // if(with_tracklets_)
    // {
    // property_map<traxel_arc_id, HypothesesGraph::base_graph>::type& traxel_arc_id_map = tracklet_graph_.get(traxel_arc_id());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = g.get(node_tracklet());
    // }


    //write this map to label file
    map<int,label_type> cplexid_label_map;
    map<int,stringstream> cplexid_label_info_map;
    map<size_t,std::vector<std::pair<size_t,double> > > cplexid_weight_class_map;

    // fill labels of Variables
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        //appearance
        for (size_t state = 0;
             state < get_model().numberOfLabels(get_appearance_node_map()[n]);
             ++state) {
            int id = clpex_variable_id_map_[make_pair(get_appearance_node_map()[n],state)];
            cplexid_label_map[id] = ((g.get(appearance_label())[n]==state)?1:0);
            LOG(logDEBUG4) <<"app\t"<< cplexid_label_map[id] << "  " << id << "  "
                           <<  get_appearance_node_map()[n] << "  " << state;
            // cplexid_label_info_map[id] <<  "# appearance id:" << id << " state:"
            // << g.get(appearance_label())[n] << "/" << state << "  node:" << get_appearance_node_map()[n]
            // << "traxel id:" <<  traxel_map[n].Id << "  ts:" << traxel_map[n].Timestep ;
            cplexid_weight_class_map[id].clear();
            cplexid_weight_class_map[id].push_back(std::make_pair(0,1.));
        }
        //disappearance
        for (size_t state = 0;
             state < get_model().numberOfLabels(get_disappearance_node_map()[n]);
             ++state) {
            int id = clpex_variable_id_map_[make_pair(get_disappearance_node_map()[n],state)];
            cplexid_label_map[id] = ((g.get(disappearance_label())[n]==state)?1:0);
            LOG(logDEBUG4) <<"dis\t"<< cplexid_label_map[id] << "  " << id << "  "
                           <<  get_disappearance_node_map()[n] << "  " << state;
            // cplexid_label_info_map[id] <<  "# disappearance id:" << id << " state:"
            // << g.get(disappearance_label())[n] << "/" << state << "  node:"
            // << get_disappearance_node_map()[n]  << "traxel id:" <<  traxel_map[n].Id << "  ts:"
            // << traxel_map[n].Timestep ;
            cplexid_weight_class_map[id].clear();
            cplexid_weight_class_map[id].push_back(std::make_pair(1,1.));
        }
        //division
        if(param_.with_divisions and get_division_node_map().count(n) != 0){

            for (size_t state = 0;
                 state < get_model().numberOfLabels(get_division_node_map()[n]);
                 ++state) {
                int id = clpex_variable_id_map_[make_pair(get_division_node_map()[n],state)];
                cplexid_label_map[id] = ((g.get(division_label())[n]==state)?1:0);
                LOG(logDEBUG4) <<"div\t"<< cplexid_label_map[id] << "  " << id << "  "
                               << get_division_node_map()[n] << "  " << state <<"   ";// << number_of_division_nodes_;
                // cplexid_label_info_map[id] <<  "# division id:" << id << " state:"
                // << g.get(division_label())[n] << "/" << state << "  node:" << div_node_map_[n]  << "traxel id:"
                // <<  traxel_map[n].Id << "  ts:" << traxel_map[n].Timestep ;
                cplexid_weight_class_map[id].clear();
                cplexid_weight_class_map[id].push_back(std::make_pair(2,1.));
            }
        }
    }
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        //move
        for (size_t state = 0; state < get_model().numberOfLabels(get_arc_map()[a]); ++state) {
            int id = clpex_variable_id_map_[make_pair(get_arc_map()[a],state)];
            cplexid_label_map[id] = ((g.get(arc_label())[a]==state)?1:0);
            LOG(logDEBUG4) <<"arc\t"<< cplexid_label_map[id] << "  " <<id << "  " <<  get_arc_map()[a] << "  " << state;
            // cplexid_label_info_map[id] <<  "# move id:" << id << " state:" << g.get(arc_label())[a] << "/"
            // << state << "  node:" << get_arc_map()[a]  << "traxel id:"
            // <<  traxel_map[g.source(a)].Id << "-->" << traxel_map[g.target(a)].Id << "  ts:"
            // << traxel_map[g.target(a)].Timestep ;

            cplexid_weight_class_map[id].clear();
            if (param_.with_tracklets and (tracklet_map[g.source(a)]).size() > 1)
            {
                cplexid_weight_class_map[id].push_back(std::make_pair(3,(tracklet_map[g.source(a)]).size()-1));
                cplexid_weight_class_map[id].push_back(std::make_pair(4,(tracklet_map[g.source(a)]).size()));
            }
            else
            {
                if (param_.with_tracklets) {
                    assert((tracklet_map[g.source(a)]).size() == 1);
                }
                cplexid_weight_class_map[id].push_back(std::make_pair(3,1.));
            }
        }
    }

    // fill labels of Factors (only second order factors need to be exported (others accounted for in variable states))
    std::map<HypothesesGraph::Node, size_t>& detection_f_node_map = get_detection_factor_node_map();

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        //detection factor detection_node_map_
        for (size_t s1 = 0; s1 < get_model().numberOfLabels(detection_f_node_map[n]); ++s1) {
            for (size_t s2 = 0; s2 < get_model().numberOfLabels(detection_f_node_map[n]); ++s2) {
                int id = clpex_factor_id_map_[make_pair(detection_f_node_map[n],make_pair(s1,s2))];
                cplexid_label_map[id] = ((g.get(appearance_label())[n]==s1 and g.get(disappearance_label())[n]==s2)?1:0);
                LOG(logDEBUG4) <<"detection\t"<< cplexid_label_map[id] <<"  "<<id<< "  "
                               <<  detection_f_node_map[n] << "  " << s1 <<"  " << s2 << endl;
                cplexid_weight_class_map[id].clear();
                if (param_.with_tracklets and (tracklet_map[n]).size() > 1)
                {
                    cplexid_weight_class_map[id].push_back(std::make_pair(3,(tracklet_map[n]).size()-1));
                    cplexid_weight_class_map[id].push_back(std::make_pair(4,(tracklet_map[n]).size()));
                }
                else{
                    if (param_.with_tracklets) {
                        assert((tracklet_map[n]).size() == 1);
                    }
                    cplexid_weight_class_map[id].push_back(std::make_pair(4,1.));
                }
                // cplexid_label_info_map[id] <<  "# factor id:" << id << " state:" << g.get(appearance_label())[n]
                // <<","<<g.get(disappearance_label())[n]<< "/" << s1 << "," << s2 << "  node:"
                // << get_appearance_node_map()[n]  << "traxel id:" <<  traxel_map[n].Id << "  ts:"
                // << traxel_map[n].Timestep ;
            }
        }
    }



    //write labels to file
    std::ofstream ground_truth_file;

    ground_truth_file.open (ground_truth_filename,std::ios::app);

    for(std::map<int,size_t>::iterator iterator = cplexid_label_map.begin(); iterator != cplexid_label_map.end(); ++iterator)
    {
        ground_truth_file << iterator->second <<"\t\t";
        for (std::vector<std::pair<size_t,double>>::iterator class_weight_pair = cplexid_weight_class_map[iterator->first].begin();
             class_weight_pair != cplexid_weight_class_map[iterator->first].end();
             ++class_weight_pair)
        {
            ground_truth_file << "#c"<< class_weight_pair->first << ":" << class_weight_pair->second << " ";
        }
        ground_truth_file<<"\tc#\tcplexid:" << iterator->first
        // << cplexid_label_info_map[iterator->first].str()
        << endl;
    }
    ground_truth_file.close();

}

} // namespace pgmlink

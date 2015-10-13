#ifndef OPENGM_UNSIGNED_INTEGER_POW_HXX_
#define OPENGM_UNSIGNED_INTEGER_POW_HXX_
#endif
#include "pgmlink/inferencemodel/constrackingexplicitinferencemodel.h"
#include <boost/python.hpp>
#include <opengm/functions/explicit_function.hxx>

namespace pgmlink
{

ConsTrackingExplicitInferenceModel::ConsTrackingExplicitInferenceModel(const Parameter& param,
        double ep_gap,
        double cplex_timeout):
    InferenceModel(param),
    number_of_transition_nodes_(0),
    number_of_division_nodes_(0),
    number_of_appearance_nodes_(0),
    number_of_disappearance_nodes_(0),
    ground_truth_filename_(""),
    inferenceWeights_(5)
{
    cplex2_param_.verbose_ = true;
    cplex2_param_.integerConstraintNodeVar_ = true;
    cplex2_param_.epGap_ = ep_gap;
    cplex2_param_.timeLimit_ = cplex_timeout;
}

void ConsTrackingExplicitInferenceModel::build_from_graph(const HypothesesGraph& hypotheses)
{

    LOG(logDEBUG) << "ConsTrackingExplicitInferenceModel::build_from_graph: add_transition_nodes";
    add_transition_nodes(hypotheses);
    LOG(logDEBUG) << "ConsTrackingExplicitInferenceModel::build_from_graph: add_appearance_nodes";
    add_appearance_nodes(hypotheses);
    LOG(logDEBUG) << "ConsTrackingExplicitInferenceModel::build_from_graph: add_disappearance_nodes";
    add_disappearance_nodes(hypotheses);

    if (param_.with_divisions)
    {
        LOG(logDEBUG) << "ConsTrackingExplicitInferenceModel::build_from_graph: add_division_nodes";
        add_division_nodes(hypotheses);
    }

    LOG(logINFO) << "number_of_transition_nodes_ = " << number_of_transition_nodes_;
    LOG(logINFO) << "number_of_appearance_nodes_ = " << number_of_appearance_nodes_;
    LOG(logINFO) << "number_of_disappearance_nodes_ = " << number_of_disappearance_nodes_;
    LOG(logINFO) << "number_of_division_nodes_ = " << number_of_division_nodes_;

    add_finite_factors(hypotheses);
    add_constraints_to_pool(hypotheses);
}

void ConsTrackingExplicitInferenceModel::fixFirstDisappearanceNodesToLabels(
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
//                constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::FixNodeValueConstraint(app_node_map_[n], appearance_labels[n]));
//                constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::FixNodeValueConstraint(dis_node_map_[n], appearance_labels[n]));
                linear_constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::FixNodeValueLinearConstraint(app_node_map_[n], appearance_labels[n]));
                linear_constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::FixNodeValueLinearConstraint(dis_node_map_[n], appearance_labels[n]));
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
//                constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::FixNodeValueConstraint(app_node_map_[n], appearance_labels[orig_n]));
//                constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::FixNodeValueConstraint(dis_node_map_[n], appearance_labels[orig_n]));
                linear_constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::FixNodeValueLinearConstraint(app_node_map_[n], appearance_labels[orig_n]));
                linear_constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::FixNodeValueLinearConstraint(dis_node_map_[n], appearance_labels[orig_n]));
            }
        }
    }
}

ConsTrackingExplicitInferenceModel::GraphicalModelType& ConsTrackingExplicitInferenceModel::get_model()
{
    return model_;
}

ConsTrackingExplicitInferenceModel::HypothesesGraphNodeMap &ConsTrackingExplicitInferenceModel::get_division_node_map()
{
    return div_node_map_;
}

ConsTrackingExplicitInferenceModel::HypothesesGraphNodeMap &ConsTrackingExplicitInferenceModel::get_appearance_node_map()
{
    return app_node_map_;
}

ConsTrackingExplicitInferenceModel::HypothesesGraphNodeMap &ConsTrackingExplicitInferenceModel::get_disappearance_node_map()
{
    return dis_node_map_;
}

ConsTrackingExplicitInferenceModel::HypothesesGraphArcMap &ConsTrackingExplicitInferenceModel::get_arc_map()
{
    return arc_map_;
}

std::map<HypothesesGraph::Node, size_t>& ConsTrackingExplicitInferenceModel::get_detection_factor_node_map()
{
    return detection_f_node_map_;
}

void ConsTrackingExplicitInferenceModel::add_appearance_nodes(const HypothesesGraph& g)
{
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        model_.addVariable(param_.max_number_objects + 1);
        app_node_map_[n] = model_.numberOfVariables() - 1;

        HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
        nodes_per_timestep_[timestep_map[n]].push_back(model_.numberOfVariables() - 1);

        assert(model_.numberOfLabels(app_node_map_[n]) == param_.max_number_objects + 1);

        property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
        ++count;
    }
    number_of_appearance_nodes_ = count;
}

void ConsTrackingExplicitInferenceModel::add_disappearance_nodes(const HypothesesGraph& g)
{
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        model_.addVariable(param_.max_number_objects + 1);
        dis_node_map_[n] = model_.numberOfVariables() - 1;

        HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
        nodes_per_timestep_[timestep_map[n]].push_back(model_.numberOfVariables() - 1);

        assert(model_.numberOfLabels(dis_node_map_[n]) == param_.max_number_objects + 1);

        property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
        ++count;
    }
    number_of_disappearance_nodes_ = count;
}

void ConsTrackingExplicitInferenceModel::add_transition_nodes(const HypothesesGraph& g)
{
    size_t count = 0;
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
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

        property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());

        ++count;
    }
    number_of_transition_nodes_ = count;
}

void ConsTrackingExplicitInferenceModel::add_division_nodes(const HypothesesGraph& g)
{
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        size_t number_of_outarcs = 0;
        for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a)
        {
            ++number_of_outarcs;
        }
        if (number_of_outarcs > 1)
        {
            model_.addVariable(2);
            div_node_map_[n] = model_.numberOfVariables() - 1;
            HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
            nodes_per_timestep_[timestep_map[n]].push_back(model_.numberOfVariables() - 1);

            assert(model_.numberOfLabels(div_node_map_[n]) == 2);

            property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
            ++count;
        }
    }
    number_of_division_nodes_ = count;
}

unsigned int ConsTrackingExplicitInferenceModel::get_number_of_division_nodes(){
    return number_of_division_nodes_;
}

void ConsTrackingExplicitInferenceModel::printResults(const HypothesesGraph& g)
{
    //very verbose print of solution
    property_map<arc_active_count, HypothesesGraph::base_graph>::type& active_arcs_count =
        g.get(arc_active_count());
    property_map<node_active_count, HypothesesGraph::base_graph>::type& active_nodes_count =
        g.get(node_active_count());
    property_map<division_active_count, HypothesesGraph::base_graph>::type& active_divisions_count =
        g.get(division_active_count());

    int c = 0;
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        c = 0;
        for( std::vector<bool>::const_iterator i = active_arcs_count[a].begin();
                i != active_arcs_count[a].end();
                ++i)
        {
            LOG(logDEBUG4) << *i;
            if (*i)
            {
                c++;
            }
        }
        LOG(logDEBUG4) << "total= " << c;
    }

    for(std::map<HypothesesGraph::Node, size_t>::const_iterator it = app_node_map_.begin();
            it != app_node_map_.end(); ++it)
    {
        c = 0;
        for( std::vector<long unsigned int>::const_iterator i = active_nodes_count[it->first].begin();
                i != active_nodes_count[it->first].end();
                ++i)
        {
            LOG(logINFO) << *i;
            if (*i > 0)
            {
                c++;
            }
        }
        LOG(logINFO) << "total= " << c << std::endl;
    }
    LOG(logDEBUG4) << "division nodes " << c << std::endl;
    for(std::map<HypothesesGraph::Node, size_t>::const_iterator it = app_node_map_.begin();
            it != app_node_map_.end(); ++it)
    {
        c = 0;
        for( std::vector<bool>::const_iterator i = active_divisions_count[it->first].begin();
                i != active_divisions_count[it->first].end();
                ++i)
        {
            LOG(logDEBUG4) << *i << " ";
            if (*i > 0)
            {
                c++;
            }
        }
        LOG(logDEBUG4) << "total= " << c << std::endl;
    }
}

size_t ConsTrackingExplicitInferenceModel::add_detection_factors(const HypothesesGraph& g, size_t factorIndex)
{

    return 0;
}

size_t ConsTrackingExplicitInferenceModel::add_transition_factors(const HypothesesGraph& g, size_t factorIndex)
{
    return 0;
}

size_t ConsTrackingExplicitInferenceModel::add_division_factors(const HypothesesGraph& g, size_t factorIndex)
{
    if(!param_.with_divisions)
    {
        return factorIndex;
    }

    return 0;
}

void ConsTrackingExplicitInferenceModel::add_finite_factors(const HypothesesGraph& g)
{
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

    LOG(logDEBUG) << "ConsTrackingExplicitInferenceModel::add_finite_factors";
    size_t factorIndex = 0;
    factorIndex = add_detection_factors(g, factorIndex);
    factorIndex = add_transition_factors(g, factorIndex);
    factorIndex = add_division_factors(g, factorIndex);
}

//set up optimizer from constraints by reading from formulated gm
void ConsTrackingExplicitInferenceModel::add_constraints_to_pool(const HypothesesGraph& g)
{
    LOG(logDEBUG) << "ConsTrackingExplicitInferenceModel::add_constraints";

//    constraint_pool_ = pgm::ConstraintPoolExplicit(param_.forbidden_cost,
//                                           param_.with_divisions,
//                                           param_.with_appearance,
//                                           param_.with_disappearance,
//                                           param_.with_misdetections_allowed);
    linear_constraint_pool_ = pgm::ConstraintPoolExplicit(param_.forbidden_cost,
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

//            constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::OutgoingConstraint(appearance_node,division_node,transition_nodes));
            linear_constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::OutgoingLinearConstraint(appearance_node,division_node,transition_nodes));
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

//            constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::IncomingConstraint(transition_nodes,disappearance_node));
            linear_constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::IncomingLinearConstraint(transition_nodes,disappearance_node));
        }

        ////
        //// disappearance/appearance coupling
        ////
        if (app_node_map_.count(n) > 0 && dis_node_map_.count(n) > 0)
        {
//            constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::DetectionConstraint((size_t)dis_node_map_[n],(size_t)app_node_map_[n]));
            linear_constraint_pool_.add_constraint(pgm::ConstraintPoolExplicit::DetectionLinearConstraint((size_t)dis_node_map_[n],(size_t)app_node_map_[n]));
        }
    }

//    constraint_pool_.force_softconstraint(!param_.with_constraints);
    //linear_constraint_pool_.force_softconstraint(!param_.with_constraints);
}

void ConsTrackingExplicitInferenceModel::set_inference_params(size_t numberOfSolutions,
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

ConsTrackingExplicitInferenceModel::IlpSolution ConsTrackingExplicitInferenceModel::infer()
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

ConsTrackingExplicitInferenceModel::IlpSolution ConsTrackingExplicitInferenceModel::extractSolution(size_t k,
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

void ConsTrackingExplicitInferenceModel::write_labeledgraph_to_file(const HypothesesGraph &g,
        const std::string &ground_truth_filename)
{
    using namespace std;

    cout << "write_labeledgraph_to_file" << endl;
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = g.get(node_tracklet());

    map<int, label_type> cplexid_label_map;
    map<int, stringstream> cplexid_label_info_map;
    map<size_t, std::vector<std::pair<size_t, double> > > cplexid_weight_class_map;

    // fill labels of Variables
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        //appearance
        for (size_t state = 0;
                state < get_model().numberOfLabels(get_appearance_node_map()[n]);
                ++state)
        {
            int id = cplex_variable_id_map_[make_pair(get_appearance_node_map()[n], state)];
            cplexid_label_map[id] = ((g.get(appearance_label())[n] == state) ? 1 : 0);
            LOG(logDEBUG4) << "app\t" << cplexid_label_map[id] << "  " << id << "  "
                           <<  get_appearance_node_map()[n] << "  " << state;
            cplexid_weight_class_map[id].clear();
            cplexid_weight_class_map[id].push_back(std::make_pair(0, 1.));
        }
        //disappearance
        for (size_t state = 0;
                state < get_model().numberOfLabels(get_disappearance_node_map()[n]);
                ++state)
        {
            int id = cplex_variable_id_map_[make_pair(get_disappearance_node_map()[n], state)];
            cplexid_label_map[id] = ((g.get(disappearance_label())[n] == state) ? 1 : 0);
            LOG(logDEBUG4) << "dis\t" << cplexid_label_map[id] << "  " << id << "  "
                           <<  get_disappearance_node_map()[n] << "  " << state;
            cplexid_weight_class_map[id].clear();
            cplexid_weight_class_map[id].push_back(std::make_pair(1, 1.));
        }
        //division
        if(param_.with_divisions and get_division_node_map().count(n) != 0)
        {

            for (size_t state = 0;
                    state < get_model().numberOfLabels(get_division_node_map()[n]);
                    ++state)
            {
                int id = cplex_variable_id_map_[make_pair(get_division_node_map()[n], state)];
                cplexid_label_map[id] = ((g.get(division_label())[n] == state) ? 1 : 0);
                LOG(logDEBUG4) << "div\t" << cplexid_label_map[id] << "  " << id << "  "
                               << get_division_node_map()[n] << "  " << state << "   "; // << number_of_division_nodes_;
                cplexid_weight_class_map[id].clear();
                cplexid_weight_class_map[id].push_back(std::make_pair(2, 1.));
            }
        }
    }
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        //move
        for (size_t state = 0; state < get_model().numberOfLabels(get_arc_map()[a]); ++state)
        {
            int id = cplex_variable_id_map_[make_pair(get_arc_map()[a], state)];
            cplexid_label_map[id] = ((g.get(arc_label())[a] == state) ? 1 : 0);
            LOG(logDEBUG4) << "arc\t" << cplexid_label_map[id] << "  " << id << "  " <<  get_arc_map()[a] << "  " << state;
            cplexid_weight_class_map[id].clear();
            if (param_.with_tracklets and (tracklet_map[g.source(a)]).size() > 1)
            {
                cplexid_weight_class_map[id].push_back(std::make_pair(3, (tracklet_map[g.source(a)]).size() - 1));
                cplexid_weight_class_map[id].push_back(std::make_pair(4, (tracklet_map[g.source(a)]).size()));
            }
            else
            {
                cplexid_weight_class_map[id].push_back(std::make_pair(3, 1.));
            }
        }
    }

    // fill labels of Factors (only second order factors need to be exported (others accounted for in variable states))
    std::map<HypothesesGraph::Node, size_t>& detection_f_node_map = get_detection_factor_node_map();

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        //detection factor detection_node_map_
        for (size_t s1 = 0; s1 < get_model().numberOfLabels(detection_f_node_map[n]); ++s1)
        {
            for (size_t s2 = 0; s2 < get_model().numberOfLabels(detection_f_node_map[n]); ++s2)
            {
                int id = cplex_factor_id_map_[make_pair(detection_f_node_map[n], make_pair(s1, s2))];
                cplexid_label_map[id] = ((g.get(appearance_label())[n] == s1 and g.get(disappearance_label())[n] == s2) ? 1 : 0);
                LOG(logDEBUG4) << "detection\t" << cplexid_label_map[id] << "  " << id << "  "
                               <<  detection_f_node_map[n] << "  " << s1 << "  " << s2 << endl;
                cplexid_weight_class_map[id].clear();
                if (param_.with_tracklets and (tracklet_map[n]).size() > 1)
                {
                    cplexid_weight_class_map[id].push_back(std::make_pair(3, (tracklet_map[n]).size() - 1));
                    cplexid_weight_class_map[id].push_back(std::make_pair(4, (tracklet_map[n]).size()));
                }
                else
                {
                    cplexid_weight_class_map[id].push_back(std::make_pair(4, 1.));
                }
            }
        }
    }

    //write labels to file
    std::ofstream ground_truth_file;

    ground_truth_file.open (ground_truth_filename, std::ios::app);

    for(std::map<int, size_t>::iterator iterator = cplexid_label_map.begin(); iterator != cplexid_label_map.end(); ++iterator)
    {
        ground_truth_file << iterator->second << "\t\t";
        for (std::vector<std::pair<size_t, double>>::iterator class_weight_pair = cplexid_weight_class_map[iterator->first].begin();
                class_weight_pair != cplexid_weight_class_map[iterator->first].end();
                ++class_weight_pair)
        {
            ground_truth_file << "#c" << class_weight_pair->first << ":" << class_weight_pair->second << " ";
        }
        ground_truth_file << "\tc#\tcplexid:" << iterator->first
                          // << cplexid_label_info_map[iterator->first].str()
                          << endl;
    }
    ground_truth_file.close();
}

void ConsTrackingExplicitInferenceModel::conclude( HypothesesGraph& g,
        HypothesesGraph& tracklet_graph,
        std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> >& tracklet2traxel_node_map,
        IlpSolution& solution)
{
    // add 'active' properties to graph
    g.add(node_active2()).add(arc_active()).add(division_active());

    property_map<node_active2, HypothesesGraph::base_graph>::type& active_nodes =
        g.get(node_active2());
    property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());
    property_map<division_active, HypothesesGraph::base_graph>::type& division_nodes =
        g.get(division_active());

    // add counting properties for analysis of perturbed models
    g.add(arc_active_count()).add(node_active_count()).add(division_active_count());

    property_map<arc_active_count, HypothesesGraph::base_graph>::type& active_arcs_count =
        g.get(arc_active_count());
    property_map<node_active_count, HypothesesGraph::base_graph>::type& active_nodes_count =
        g.get(node_active_count());
    property_map<division_active_count, HypothesesGraph::base_graph>::type& active_divisions_count =
        g.get(division_active_count());

    if (!param_.with_tracklets)
    {
        tracklet_graph.add(tracklet_intern_arc_ids()).add(traxel_arc_id());
    }
    property_map<tracklet_intern_arc_ids, HypothesesGraph::base_graph>::type& tracklet_arc_id_map =
        tracklet_graph.get(tracklet_intern_arc_ids());
    property_map<traxel_arc_id, HypothesesGraph::base_graph>::type& traxel_arc_id_map =
        tracklet_graph.get(traxel_arc_id());

    int iterStep = active_nodes_count[get_appearance_node_map().begin()->first].size();
    bool isMAP = (iterStep == 0);

    if (isMAP)
    {
        //initialize vectors for storing optimizer results
        for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
        {
            active_arcs_count.set(a, std::vector<bool>());
        }
        for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
        {
            active_nodes_count.set(n, std::vector<long unsigned int>());
            active_divisions_count.set(n, std::vector<bool>());
        }
    }

    //initialize node counts by 0
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        active_nodes_count.get_value(n).push_back(0);
        active_divisions_count.get_value(n).push_back(0);
    }

    //initialize arc counts by 0
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        active_arcs.set(a, false);
        active_arcs_count.get_value(a).push_back(0);
    }

    // write state after inference into 'active'-property maps
    // the node is also active if its appearance node is active
    for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = get_appearance_node_map().begin();
            it != get_appearance_node_map().end();
            ++it)
    {
        if (param_.with_tracklets)
        {
            // set state of tracklet nodes
            std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map[it->first];

            for (std::vector<HypothesesGraph::Node>::const_iterator tr_n_it = traxel_nodes.begin();
                    tr_n_it != traxel_nodes.end();
                    ++tr_n_it)
            {
                HypothesesGraph::Node n = *tr_n_it;
                active_nodes.set(n, solution[it->second]);
                active_nodes_count.get_value(n)[iterStep] = solution[it->second];
            }

            // set state of tracklet internal arcs
            std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
            for (std::vector<int>::const_iterator arc_id_it = arc_ids.begin();
                    arc_id_it != arc_ids.end(); ++arc_id_it)
            {
                HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
                assert(active_arcs[a] == false);
                if (solution[it->second] > 0)
                {

                    active_arcs.set(a, true);
                    active_arcs_count.get_value(a)[iterStep] = true;

                    assert(active_nodes[g.source(a)] == solution[it->second]
                           && "tracklet internal arcs must have the same flow as their connected nodes");
                    assert(active_nodes[g.target(a)] == solution[it->second]
                           && "tracklet internal arcs must have the same flow as their connected nodes");
                }
            }
        }
        else
        {
            active_nodes.set(it->first, solution[it->second]);
            active_nodes_count.get_value(it->first)[iterStep] = solution[it->second];
        }
    }
    // the node is also active if its disappearance node is active
    for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = get_disappearance_node_map().begin();
            it != get_disappearance_node_map().end(); ++it)
    {
        if (solution[it->second] > 0)
        {
            if (param_.with_tracklets)
            {
                // set state of tracklet nodes
                std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map[it->first];
                for (std::vector<HypothesesGraph::Node>::const_iterator tr_n_it =
                            traxel_nodes.begin(); tr_n_it != traxel_nodes.end(); ++tr_n_it)
                {
                    HypothesesGraph::Node n = *tr_n_it;

                    if (active_nodes[n] == 0)
                    {
                        active_nodes.set(n, solution[it->second]);
                        active_nodes_count.get_value(n)[iterStep] = solution[it->second];

                    }
                    else
                    {
                        assert(active_nodes[n] == solution[it->second]);
                    }
                }
                // set state of tracklet internal arcs
                std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
                for (std::vector<int>::const_iterator arc_id_it = arc_ids.begin();
                        arc_id_it != arc_ids.end(); ++arc_id_it)
                {
                    HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
                    if (solution[it->second] > 0)
                    {

                        active_arcs.set(a, true);
                        active_arcs_count.get_value(a)[iterStep] = true;

                        assert(active_nodes[g.source(a)] == solution[it->second]
                               && "tracklet internal arcs must have the same flow as their connected nodes");
                        assert(active_nodes[g.target(a)] == solution[it->second]
                               && "tracklet internal arcs must have the same flow as their connected nodes");
                    }
                }
            }
            else
            {

                if (active_nodes[it->first] == 0)
                {
                    active_nodes.set(it->first, solution[it->second]);
                    active_nodes_count.get_value(it->first)[iterStep] = solution[it->second];

                }
                else
                {
                    assert(active_nodes[it->first] == solution[it->second]);
                }
            }
        }
    }

    for (std::map<HypothesesGraph::Arc, size_t>::const_iterator it = get_arc_map().begin();
            it != get_arc_map().end(); ++it)
    {
        if (solution[it->second] >= 1)
        {
            if (param_.with_tracklets)
            {
                active_arcs.set(g.arcFromId((traxel_arc_id_map[it->first])), true);
                active_arcs_count.get_value(g.arcFromId((traxel_arc_id_map[it->first])))[iterStep] = true;
            }
            else
            {
                active_arcs.set(it->first, true);
                active_arcs_count.get_value(it->first)[iterStep] = true;
            }
        }
    }
    // write division node map
    if (param_.with_divisions)
    {
        for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = get_division_node_map().begin();
                it != get_division_node_map().end(); ++it)
        {
            division_nodes.set(it->first, false);
        }
        for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = get_division_node_map().begin();
                it != get_division_node_map().end(); ++it)
        {

            if (solution[it->second] >= 1)
            {
                if (param_.with_tracklets)
                {
                    // set division property for the last node in the tracklet
                    HypothesesGraph::Node n = tracklet2traxel_node_map[it->first].back();
                    division_nodes.set(n, true);
                    active_divisions_count.get_value(n)[iterStep] = true;
                }
                else
                {
                    division_nodes.set(it->first, true);
                    active_divisions_count.get_value(it->first)[iterStep] = true;
                }
            }
        }
    }
}

ConsTrackingExplicitInferenceModel::GraphicalModelType ConsTrackingExplicitInferenceModel::model(){
    return model_;
}

} // namespace pgmlink

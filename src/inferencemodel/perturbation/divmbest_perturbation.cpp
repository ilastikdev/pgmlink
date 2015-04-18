#include "pgmlink/inferencemodel/perturbation/divmbest_perturbation.h"

namespace pgmlink
{

DivMBestPerturbation::DivMBestPerturbation(const Parameter &perturbation_param,
                           const InferenceModel::Parameter &inf_param):
    Perturbation(perturbation_param, inf_param)
{
    if(perturbation_param_.distributionId != DiverseMbest)
        throw std::runtime_error("Cannot construct DivMBestPerturbation "
                                 "when another distribution is specified");
}

void DivMBestPerturbation::push_away_from_solution(const PertGmType &model, std::vector<size_t> solution)
{
    size_t num_factors = model.numberOfFactors();
    if(deterministic_offset_.size() != num_factors)
        deterministic_offset_.resize(num_factors);

    for (size_t factorId = 0; factorId < num_factors; ++factorId)
    {
        PertGmType::FactorType factor = model[factorId];
        vector<size_t> varIndices;
        for (PertGmType::FactorType::VariablesIteratorType ind = factor.variableIndicesBegin();
                ind != factor.variableIndicesEnd();
                ++ind)
        {
            varIndices.push_back(solution[*ind]);
        }
        deterministic_offset_[factorId].push_back(varIndices);
    }
}

double DivMBestPerturbation::generateRandomOffset(EnergyType energyIndex,
                                          double energy,
                                          Traxel tr,
                                          Traxel tr2,
                                          size_t state,
                                          boost::shared_ptr<InferenceModel::TransitionPredictionsMap> transition_predictions)
{

    LOG(logDEBUG4) << "generateRandomOffset()";

    LOG(logDEBUG4) << "DiverseMBest/MBestCPLEX: random offset 0";
    return 0;
}

size_t DivMBestPerturbation::add_div_m_best_perturbation(marray::Marray<double>& energies,
        EnergyType energy_type,
        size_t factorIndex)
{
    if (perturbation_param_.distributionId == DiverseMbest)
    {
        std::vector<std::vector<size_t> >& indexlist = deterministic_offset_[factorIndex];
        for (std::vector<std::vector<size_t> >::iterator index = indexlist.begin(); index != indexlist.end(); index++)
        {
            energies(index->begin()) += perturbation_param_.distributionParam[energy_type];
        }
        factorIndex++;
    }
    return factorIndex;
}

void DivMBestPerturbation::registerOriginalGraph(const HypothesesGraph* g,
                                                 std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> >* tracklet2traxel_node_map)
{
    original_graph_ = g;
    tracklet2traxel_node_map_ = tracklet2traxel_node_map;
}

double DivMBestPerturbation::getDivMBestOffset(EnergyType energy_type,
                                               const HypothesesGraph& g,
                                               HypothesesGraph::Node n,
                                               HypothesesGraph::Arc a,
                                               size_t state)
{
    double energy = 0.0;

    if(!param_.with_tracklets)
    {
        switch(energy_type)
        {
            case Detection:
            {
                assert(g.has_property(node_active_count()));
                property_map<node_active_count, HypothesesGraph::base_graph>::type& active_nodes_count = g.get(node_active_count());

                for(size_t label : active_nodes_count[n])
                {
                    if(label == state)
                    {
                        energy += perturbation_param_.distributionParam[energy_type];
                    }
                }
            } break;
            case Division:
            {
                assert(g.has_property(division_active_count()));
                property_map<division_active_count, HypothesesGraph::base_graph>::type& active_division_count = g.get(division_active_count());

                for(bool label : active_division_count[n])
                {
                    if(label == bool(state))
                    {
                        energy += perturbation_param_.distributionParam[energy_type];
                    }
                }
            } break;
            case Transition:
            {
                // TODO: push away from actual number of cells going along this arc, as for CPLEX above?
                assert(g.has_property(arc_active_count()));
                property_map<arc_active_count, HypothesesGraph::base_graph>::type& active_arc_count = g.get(arc_active_count());

                for(bool label : active_arc_count[a])
                {
                    if(label == bool(state))
                    {
                        energy += perturbation_param_.distributionParam[energy_type];
                    }
                }
            } break;
        }
    }
    else
    {
        switch(energy_type)
        {
            case Detection:
            {
                // it should not matter which one we use, they should all have the same label
                // but the last one can also divide
                HypothesesGraph::Node orig_n = (*tracklet2traxel_node_map_)[n].back();

                assert(original_graph_->has_property(node_active_count()));
                property_map<node_active_count, HypothesesGraph::base_graph>::type& active_nodes_count = original_graph_->get(node_active_count());

                for(size_t label : active_nodes_count[orig_n])
                {
                    if(label == state)
                    {
                        energy += perturbation_param_.distributionParam[energy_type];
                    }
                }
            } break;
            case Division:
            {
                // it should not matter which one we use, they should all have the same label
                // but the last one can also divide
                HypothesesGraph::Node orig_n = (*tracklet2traxel_node_map_)[n].back();

                assert(original_graph_->has_property(division_active_count()));
                property_map<division_active_count, HypothesesGraph::base_graph>::type& active_division_count = original_graph_->get(division_active_count());

                for(bool label : active_division_count[orig_n])
                {
                    if(label == bool(state))
                    {
                        energy += perturbation_param_.distributionParam[energy_type];
                    }
                }
            } break;
            case Transition:
            {
                // get original arc via traxel_arc_id
                HypothesesGraph::Arc orig_a = original_graph_->arcFromId(g.get(traxel_arc_id())[a]);

                // TODO: push away from actual number of cells going along this arc, as for CPLEX above?
                assert(original_graph_->has_property(arc_active_count()));
                property_map<arc_active_count, HypothesesGraph::base_graph>::type& active_arc_count = original_graph_->get(arc_active_count());

                for(bool label : active_arc_count[orig_a])
                {
                    if(label == bool(state))
                    {
                        energy += perturbation_param_.distributionParam[energy_type];
                    }
                }
            } break;
        }
    }

    return energy;
}

} // namespace pgmlink

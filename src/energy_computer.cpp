#include "pgmlink/energy_computer.h"
#include "pgmlink/features/featurestore.h"
#include <limits>

namespace pgmlink
{

EnergyComputer::EnergyComputer(
	Parameter& param,
	bool convexify,
	const std::string& detectionEnergyName,
	const std::string& divisionEnergyName,
	const std::string& appearanceEnergyName,
	const std::string& disappearanceEnergyName,
	const std::string& transitionEnergyName):
	param_(param),
	convexify_(convexify)
{
	energyNames_[Appearance] = appearanceEnergyName;
	energyNames_[Disappearance] = disappearanceEnergyName;
	energyNames_[Division] = divisionEnergyName;
	energyNames_[Transition] = transitionEnergyName;
	energyNames_[Detection] = detectionEnergyName;
}

void EnergyComputer::operator()(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs)
{
	for(HypothesesGraph::NodeIt n(graph); n != lemon::INVALID; ++n)
	{
		computeDetectionEnergy(graph, fs, n);
		computeDivisionEnergy(graph, fs, n);
		computeAppearanceEnergy(graph, fs, n);
		computeDisappearanceEnergy(graph, fs, n);
	}

	for(HypothesesGraph::ArcIt a(graph); a != lemon::INVALID; ++a)
	{
		computeTransitionEnergy(graph, fs, a);
	}
}

void EnergyComputer::storeNodeEnergies(
	HypothesesGraph& graph, 
	boost::shared_ptr<FeatureStore> fs, 
	HypothesesGraph::Node n, 
	EnergyType t,
	const feature_array& energies)
{
	TraxelMap& traxel_map = graph.get(node_traxel());
    TrackletMap& tracklet_map = graph.get(node_tracklet());

    if(param_.with_tracklets)
    {
    	for(std::vector<Traxel>::const_iterator trax_it = tracklet_map[n].begin();
            trax_it != tracklet_map[n].end(); 
            ++trax_it)
		{
			const Traxel& tr = *trax_it;
		    assert(tr.get_feature_store() == fs);
		    fs->get_traxel_features(tr)[energyNames_[t]] = energies;
		}
    }
    else
    {
    	const Traxel& tr = traxel_map[n];
	    assert(tr.get_feature_store() == fs);
	    fs->get_traxel_features(tr)[energyNames_[t]] = energies;
	}
}

void EnergyComputer::convexifyEnergies(feature_array& energies)
{
	if(convexify_)
	{
		// TODO: implement me!
	}
}

void EnergyComputer::computeDetectionEnergy(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs, HypothesesGraph::Node n)
{
	TraxelMap& traxel_map = graph.get(node_traxel());
    TrackletMap& tracklet_map = graph.get(node_tracklet());

    feature_array energyPerCellCount(param_.max_number_objects, std::numeric_limits<feature_type>::infinity());
    for(size_t state = 0; state < param_.max_number_objects; ++state)
    {
    	feature_type& energy = energyPerCellCount[state];

    	if(param_.with_tracklets)
        {
            energy = 0;
            // add all detection factors of the internal nodes
            for(std::vector<Traxel>::const_iterator trax_it = tracklet_map[n].begin();
                trax_it != tracklet_map[n].end(); 
                ++trax_it)
            {
                energy += param_.detection(*trax_it, state);
            }

            // add all transition factors of the internal arcs
            Traxel tr_prev;
            bool first = true;
            for(std::vector<Traxel>::const_iterator trax_it = tracklet_map[n].begin();
                trax_it != tracklet_map[n].end(); 
                ++trax_it)
            {
                LOG(logDEBUG4) << "internal arcs traxel " << *trax_it;
                Traxel tr = *trax_it;
                if(!first)
                    energy += param_.transition( get_transition_probability(tr_prev, tr, state) );
//                    energy += param_.transition( tr_prev, tr, state);
                else
                    first = false;
                tr_prev = tr;
            }
        }
        else
            energy = param_.detection(traxel_map[n], state);
    }

	convexifyEnergies(energyPerCellCount);
    storeNodeEnergies(graph, fs, n, Detection, energyPerCellCount);
}

void EnergyComputer::computeDivisionEnergy(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs, HypothesesGraph::Node n)
{
	TraxelMap& traxel_map = graph.get(node_traxel());
    TrackletMap& tracklet_map = graph.get(node_tracklet());

	Traxel tr;
    if (param_.with_tracklets)
    {
        tr = tracklet_map[n].back();
    }
    else
    {
        tr = traxel_map[n];
    }

    feature_array division_energy;
    for (size_t state = 0; state <= 1; ++state)
    {
        division_energy.push_back(param_.division(tr, state));
    }

	// no need to convexify as this is binary
    storeNodeEnergies(graph, fs, n, Division, division_energy);
}

void EnergyComputer::computeAppearanceEnergy(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs, HypothesesGraph::Node n)
{
	TraxelMap& traxel_map = graph.get(node_traxel());
    TrackletMap& tracklet_map = graph.get(node_tracklet());

	int node_begin_time = -1;
    if (param_.with_tracklets)
    {
        node_begin_time = tracklet_map[n].front().Timestep;
    }
    else
    {
        node_begin_time = traxel_map[n].Timestep;
    }

    double energy;
    // "<" holds if there are only tracklets in the first frame
    if (node_begin_time <= graph.earliest_timestep())
    {
        energy = 0.0;
    }
    else
    {
        if (param_.with_tracklets)
        {
            energy = param_.appearance_cost_fn(tracklet_map[n].front());
        }
        else
        {
            energy = param_.appearance_cost_fn(traxel_map[n]);
        }
    }

    feature_array energyPerCellCount;
    for(size_t state = 0; state < param_.max_number_objects; ++state)
    	energyPerCellCount.push_back(energy * state);

	convexifyEnergies(energyPerCellCount);
    storeNodeEnergies(graph, fs, n, Appearance, energyPerCellCount);
}

void EnergyComputer::computeDisappearanceEnergy(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs, HypothesesGraph::Node n)
{
	TraxelMap& traxel_map = graph.get(node_traxel());
    TrackletMap& tracklet_map = graph.get(node_tracklet());

    int node_end_time = -1;
    if (param_.with_tracklets)
    {
        node_end_time = tracklet_map[n].back().Timestep;
    }
    else
    {
        node_end_time = traxel_map[n].Timestep;
    }

    double energy;
    // "<" holds if there are only tracklets in the last frame
    if (node_end_time < graph.latest_timestep())
    {
        if (param_.with_tracklets)
        {
            energy = param_.disappearance_cost_fn(tracklet_map[n].back());
        }
        else
        {
            energy = param_.disappearance_cost_fn(traxel_map[n]);
        }
    }
    else
    {
        energy = 0.0;
    }

    feature_array energyPerCellCount;
    for(size_t state = 0; state < param_.max_number_objects; ++state)
    	energyPerCellCount.push_back(energy * state);

	convexifyEnergies(energyPerCellCount);
    storeNodeEnergies(graph, fs, n, Disappearance, energyPerCellCount);
}

double EnergyComputer::get_transition_prob(double distance, size_t state, double alpha) const
{
    double prob = exp(-distance / alpha);
    if (state == 0)
    {
        return 1 - prob;
    }
    return prob;
}

double EnergyComputer::get_transition_probability(Traxel& tr1, Traxel& tr2, size_t state) const
{
    double prob;
    double distance = 0;
    if (param_.with_optical_correction)
    {
        distance = tr1.distance_to_corr(tr2);
    }
    else
    {
        distance = tr1.distance_to(tr2);
    }
    prob = get_transition_prob(distance, state, param_.transition_parameter);
    LOG(logDEBUG4) << "get_transition_probability(): using deterministic function: " << tr1
                   << " " << tr2 << " [" << state << "] = " << prob << "; distance = " << distance;
    assert(prob >= 0 && prob <= 1);
    return prob;
}

void EnergyComputer::computeTransitionEnergy(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs, HypothesesGraph::Arc a)
{
	TraxelMap& traxel_map = graph.get(node_traxel());
    TrackletMap& tracklet_map = graph.get(node_tracklet());

    Traxel tr1, tr2;
    if (param_.with_tracklets)
    {
        tr1 = tracklet_map[graph.source(a)].back();
        tr2 = tracklet_map[graph.target(a)].front();
    }
    else
    {
        tr1 = traxel_map[graph.source(a)];
        tr2 = traxel_map[graph.target(a)];
    }

    feature_array energyPerCellCount;
    for (size_t state = 0; state <= param_.max_number_objects; ++state)
    {
        energyPerCellCount.push_back(param_.transition(get_transition_probability(tr1, tr2, state)));
//        energyPerCellCount.push_back(param_.transition(tr1, tr2, state));
    }

	convexifyEnergies(energyPerCellCount);
	fs->get_traxel_features(tr1, tr2)[energyNames_[Transition]] = energyPerCellCount;
}
	
} // end namespace pgmlink

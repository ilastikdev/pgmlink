#ifndef ENERGY_COMPUTER_H
#define ENERGY_COMPUTER_H

#include <iostream>
#include <map>
#include <boost/shared_ptr.hpp>
#include <opengm/graphicalmodel/graphicalmodel.hxx>

#include "conservationtracking_parameter.h"
#include "hypotheses.h"
#include "inferencemodel/inferencemodel.h" // for EnergyType enum

namespace pgmlink
{

class FeatureStore;

/// This class takes a hypotheses graph, a parameter object and a feature store,
/// and evaluates the energy/cost functions for each detection,appearance,disappearance,division and transition.
/// It can be queried for the name of each respective event's cost (which references it in the FeatureStore),
/// and also returns the weight of each of these. 
/// An InferenceModel then just needs to get the event's cost (TODO: later it also needs to multiply it by the weight).
/// The energy can also be convexified in this step.
///
/// ATTENTION: it does NOT add noise to the energies for perturbations!
class EnergyComputer
{
public:
	/// set up and configure
	EnergyComputer(
		Parameter& param,
		bool convexify = false,
		const std::string& detectionEnergyName = "detEnergy",
		const std::string& divisionEnergyName = "divEnergy",
		const std::string& appearanceEnergyName = "appEnergy",
		const std::string& disappearanceEnergyName = "disEnergy",
		const std::string& transitionEnergyName = "transEnergy");

	/// perform the computation and store results in featurestore
	void operator()(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs);

	/// get the feature name used in the FeatureStore so that it can be queried lateron
	const std::string getEnergyName(EnergyType t) const { return energyNames_.at(t); }
	// const double getEnergyWeight(t) const { return energyWeights_.at(t); }

private:
	typedef property_map<node_traxel, HypothesesGraph::base_graph>::type TraxelMap;
	typedef property_map<node_tracklet, HypothesesGraph::base_graph>::type TrackletMap;

private:
	/// Compute and store the energy of an event in the respective traxel.
	/// In the case of tracklets, the energy is stored in all contained traxels for simplicity.
	void computeDetectionEnergy(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs, HypothesesGraph::Node n);
	void computeDivisionEnergy(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs, HypothesesGraph::Node n);
	void computeAppearanceEnergy(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs, HypothesesGraph::Node n);
	void computeDisappearanceEnergy(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs, HypothesesGraph::Node n);
	void computeTransitionEnergy(HypothesesGraph& graph, boost::shared_ptr<FeatureStore> fs, HypothesesGraph::Arc a);
	void storeNodeEnergies(
		HypothesesGraph& graph, 
		boost::shared_ptr<FeatureStore> fs, 
		HypothesesGraph::Node n, 
		EnergyType t,
		const feature_array& energies);
	void convexifyEnergies(feature_array& energies);
	double get_transition_probability(Traxel& tr1, Traxel& tr2, size_t state) const;
	double get_transition_prob(double distance, size_t state, double alpha) const;

private:
	Parameter& param_;
	bool convexify_;
	
	std::map<EnergyType, std::string> energyNames_;
	// std::map<EnergyType, double> energyWeights_;
};

} // end namespace pgmlink

#endif // ENERGY_COMPUTER_H

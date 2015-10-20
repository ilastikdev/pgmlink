#ifndef DIV_M_BEST_PERTURBATION_H
#define DIV_M_BEST_PERTURBATION_H

//Random distribs
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/math/distributions/normal.hpp>

#include <vector>
#include <iostream>

#include <boost/function.hpp>

#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/inference.hxx>

#include "pgmlink/traxels.h"
#include "pgmlink/uncertaintyParameter.h"
#include "pgmlink/inferencemodel/inferencemodel.h"
#include "pgmlink/inferencemodel/perturbation/perturbation.h"

namespace pgmlink
{

class DivMBestPerturbation : public Perturbation
{
public: // API
    DivMBestPerturbation(const Perturbation::Parameter& perturbation_param, const pgmlink::Parameter& inf_param);


    virtual void push_away_from_solution(const PertGmType& model, std::vector<size_t> solution);

    virtual double generateRandomOffset(EnergyType parameterIndex,
                                        double energy = 0,
                                        Traxel tr = 0,
                                        Traxel tr2 = 0,
                                        size_t state = 0,
                                        boost::shared_ptr<InferenceModel::TransitionPredictionsMap> transition_predictions
                                            = boost::shared_ptr<InferenceModel::TransitionPredictionsMap>());

    virtual size_t add_div_m_best_perturbation(marray::Marray<double> &energies,
                                               EnergyType energy_type,
                                               size_t factorIndex);

    virtual double getDivMBestOffset(EnergyType energy_type,
                                     const HypothesesGraph &g,
                                     HypothesesGraph::Node n,
                                     HypothesesGraph::Arc a,
                                     size_t state);

    void registerOriginalGraph(const HypothesesGraph *g,
                               std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > *tracklet2traxel_node_map);
protected: // Methods
    double weightedNegLog(double prob, double weight);
    double inverseWeightedNegLog(double energy, double weight);
    double sigmoid(double x);
    double inverse_sigmoid(double x);
    double get_transition_variance(Traxel &tr1,
                                   Traxel &tr2,
                                   boost::shared_ptr<InferenceModel::TransitionPredictionsMap> transition_predictions);
    double sample_with_classifier_variance(double mean, double variance);

protected: // Members
    DeterministicOffset deterministic_offset_;
    const HypothesesGraph* original_graph_;
    std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> >* tracklet2traxel_node_map_;
};

} // namespace pgmlink

#endif // PERTURBATION_H

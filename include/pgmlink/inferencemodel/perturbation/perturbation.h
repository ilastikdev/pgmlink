#ifndef PERTURBATION_H
#define PERTURBATION_H

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

namespace pgmlink
{

class Perturbation
{
public: // Typedefs
    typedef std::vector< std::vector< std::vector<size_t> > > DeterministicOffset;
    typedef typename boost::variate_generator<boost::mt19937, boost::normal_distribution<> > normalRNGType;
    typedef typename boost::variate_generator<boost::mt19937, boost::uniform_real<> > uniformRNGType;

public: // Parameter object
    class Parameter
    {
    public:
        DistrId distributionId;
        std::vector<double> distributionParam;

        double division_weight;
        double detection_weight;
        double transition_weight;
    };

public: // API
    Perturbation(const Parameter& perturbation_param, const pgmlink::Parameter& inf_param);

    virtual double generateRandomOffset(EnergyType parameterIndex,
                                        double energy = 0,
                                        Traxel tr = 0,
                                        Traxel tr2 = 0,
                                        size_t state = 0,
                                        boost::shared_ptr<InferenceModel::TransitionPredictionsMap> transition_predictions
                                            = boost::shared_ptr<InferenceModel::TransitionPredictionsMap>()) = 0;

    virtual size_t add_div_m_best_perturbation(marray::Marray<double> &energies,
                                               EnergyType energy_type,
                                               size_t factorIndex);

    virtual double getDivMBestOffset(EnergyType energy_type,
                                     const HypothesesGraph& g,
                                     HypothesesGraph::Node n,
                                     HypothesesGraph::Arc a,
                                     size_t state);

protected: // Members
    boost::math::normal standard_gaussian_distribution;
    Parameter perturbation_param_;
    pgmlink::Parameter param_;

    boost::mt19937 rng_;
    normalRNGType random_normal_;
    uniformRNGType random_uniform_;
};

} // namespace pgmlink

#endif // PERTURBATION_H

#ifndef PERTURBEDINFERENCEMODEL_H
#define PERTURBEDINFERENCEMODEL_H

//Random distribs
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/math/distributions/normal.hpp>


#include <vector>
#include <iostream>
#include <opengm/datastructures/marray/marray.hxx>
#include "constrackinginferencemodel.h"
#include "uncertaintyParameter.h"

namespace pgmlink
{

/**
 * @brief The PerturbedInferenceModel class specializes the ConsTrackingInferenceModel
 * and adds perturbation to the unaries by providing non-zero perturbations in virtual methods.
 */
class PerturbedInferenceModel : public ConsTrackingInferenceModel
{
public: // Typedefs
    typedef std::vector< std::vector< std::vector<size_t> > > DeterministicOffset;
    typedef typename boost::variate_generator<boost::mt19937, boost::normal_distribution<> > normalRNGType;
    typedef typename boost::variate_generator<boost::mt19937, boost::uniform_real<> > uniformRNGType;

public: // Parameter object
    class Parameter{
    public:
        DistrId distributionId;
        std::vector<double> distributionParam;

        double division_weight;
        double detection_weight;
        double transition_weight;
    };

public: // API
    PerturbedInferenceModel(
            const ConsTrackingInferenceModel::Parameter& param,
            const Parameter& perturbation_param,
            double ep_gap,
            double cplex_timeout);

    void perturb(DeterministicOffset* det_off);

protected: // Methods
    double weightedNegLog(double prob, double weight);
    double inverseWeightedNegLog(double energy, double weight);
    double sigmoid(double x);
    double inverse_sigmoid(double x);
    double get_transition_variance(Traxel &tr1, Traxel &tr2);
    double sample_with_classifier_variance(double mean, double variance);

    // specialize virtual methods
    virtual double generateRandomOffset(EnergyType parameterIndex,  double energy = 0, Traxel tr = 0, Traxel tr2 = 0, size_t state = 0);
    virtual size_t add_div_m_best_perturbation(marray::Marray<double> &energies, EnergyType energy_type, size_t factorIndex);

protected: // Members
    boost::math::normal standard_gaussian_distribution;
    DeterministicOffset* deterministic_offset_;
    Parameter perturbation_param_;

    boost::mt19937 rng_;
    normalRNGType random_normal_;
    uniformRNGType random_uniform_;
};

} // namespace pgmlink

#endif // PERTURBEDINFERENCEMODEL_H

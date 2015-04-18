#include "pgmlink/inferencemodel/perturbation/perturbation.h"

namespace pgmlink
{

Perturbation::Perturbation(const Parameter &perturbation_param,
                           const InferenceModel::Parameter &inf_param):
    standard_gaussian_distribution(0.0, 1.0),
    rng_(42),
    random_normal_(rng_, boost::normal_distribution<>(0, 1)),
    random_uniform_(rng_, boost::uniform_real<>(0, 1)),
    perturbation_param_(perturbation_param),
    param_(inf_param)
{

}

size_t Perturbation::add_div_m_best_perturbation(marray::Marray<double>& energies,
        EnergyType energy_type,
        size_t factorIndex)
{
    return factorIndex;
}

double Perturbation::getDivMBestOffset(EnergyType energy_type,
                                       const HypothesesGraph &g,
                                       HypothesesGraph::Node n,
                                       HypothesesGraph::Arc a,
                                       size_t state)
{
    return 0.0;
}

} // namespace pgmlink

#include "pgmlink/perturbedinferencemodel.h"

namespace pgmlink
{

PerturbedInferenceModel::PerturbedInferenceModel(
        const ConsTrackingInferenceModel::Parameter& param,
        const Parameter &perturbation_param):
    ConsTrackingInferenceModel(param),
    standard_gaussian_distribution(0.0, 1.0),
    rng_(42),
    random_normal_(rng_,boost::normal_distribution<>(0, 1)),
    random_uniform_(rng_,boost::uniform_real<>(0,1)),
    deterministic_offset_(NULL),
    perturbation_param_(perturbation_param)
{
}

double PerturbedInferenceModel::get_transition_variance(Traxel& tr1, Traxel& tr2) {
    double var;

//    if (transition_classifier_.ptr()==boost::python::object().ptr()){
//        var = perturbation_param_.distributionParam[Transition];
//        LOG(logDEBUG4) << "using constant transition variance " << var;
//        if (var < 0) {
//            throw std::runtime_error("the transition variance must be positive");
//        }
//    } else {
//        TransitionPredictionsMap::const_iterator it = transition_predictions_.find(std::make_pair(tr1, tr2));
//        if ( it == transition_predictions_.end() ) {
//            throw std::runtime_error("cannot find prob/var. get_transition_probability must be called first");
//        }
//        var = it->second.second;
//        LOG(logDEBUG4) << "using GPC transition variance " << var;
//    }

    return var;
}

double PerturbedInferenceModel::sigmoid(double x){
    return cdf(standard_gaussian_distribution, x);
}

double PerturbedInferenceModel::inverse_sigmoid(double x){
    if (x <= 0.000001) {
        x = 0.000001;
    }
    assert(x > 0 && x <= 1);
    double res = 0.;
    try {
        res = quantile(standard_gaussian_distribution, x);
    } catch (...) {
        LOG(logERROR) << "error in inverse_sigmoid(), using 0.; x = " << x;
    }
    return res;
}

double PerturbedInferenceModel::sample_with_classifier_variance(double mean, double variance){
    //Map probability through inverse sigmoid to recover the classifier prediction.
    //Then use the corresponding variance stored in the traxel to sample from
    //a gaussian distribution. Map the result back through the sigmoid to obtain
    //perturbed probability, from which we finally compute the offset.
    double variance_factor = sqrt(1+variance);
    //apply inverse_sigmoid
    double mean_recovered = inverse_sigmoid(mean)*variance_factor;
    double new_sample = random_normal_()*variance+mean_recovered;
    double new_probability =  sigmoid(new_sample/variance_factor);
    return new_probability;
}


double PerturbedInferenceModel::inverseWeightedNegLog(double energy, double weight) {
    assert (weight != 0);
    return exp(-energy/weight);
}

double PerturbedInferenceModel::weightedNegLog(double prob, double weight) {
    if (prob <= 0.00001) {
        prob = 0.00001;
    }
    return -weight * log(prob);
}

void PerturbedInferenceModel::perturb(PerturbedInferenceModel::DeterministicOffset *det_off)
{
    deterministic_offset_ = det_off;
}

double PerturbedInferenceModel::generateRandomOffset(EnergyType energyIndex, double energy, Traxel tr, Traxel tr2, size_t state) {

    LOG(logDEBUG4) << "generateRandomOffset()";

    double rand;

    switch (perturbation_param_.distributionId) {
        case GaussianPertubation: //normal distribution
            LOG(logDEBUG4) << "GaussianPerturbation";
            if (energyIndex >= perturbation_param_.distributionParam.size()) {
                throw std::runtime_error("sigma is not set correctly");
            }
            return random_normal_() * perturbation_param_.distributionParam[energyIndex];
        case ClassifierUncertainty://sample from Gaussian Distribution where variance comes from Classifier
            LOG(logDEBUG4) << "ClassifierUncertainty";
            {
                double mean, variance, perturbed_mean, new_energy_offset;
                FeatureMap::const_iterator it;

                switch (energyIndex) {
                    case Detection:
                        if (perturbation_param_.detection_weight == 0) {
                            return 0.;
                        }
                        // this assumes that we use the NegLog as energy function
                        // TODO: write an inverse() function for each energy function
                        mean = inverseWeightedNegLog(energy, perturbation_param_.detection_weight); // convert energy to probability
                        it = tr.features.find("detProb_Var");
                        if (it == tr.features.end()) {
                            throw std::runtime_error("feature detProb_Var does not exist");
                        }
                        if (param_.max_number_objects == 1) {
                            // in the binary classification case we only have one variance value
                            variance = it->second[0];
                        } else {
                            // one vs. all has N variance values
                            variance = it->second[state];
                        }
                        if (variance == 0) {
                            // do not perturb
                            LOG(logDEBUG3) << "Detection: variance 0. -> do not perturb";
                            return 0.;
                        }
                        perturbed_mean = sample_with_classifier_variance(mean,variance);
                        // FIXME: the respective energy function should be used here
                        new_energy_offset = weightedNegLog(perturbed_mean, perturbation_param_.detection_weight) - energy;
                        LOG(logDEBUG3) << "Detection: old energy: " << energy << "; new energy offset: " << new_energy_offset;
                        return new_energy_offset;
                    case Division:
                        if (perturbation_param_.division_weight == 0) {
                            return 0.;
                        }
                        // this assumes that we use the NegLog as energy function
                        // TODO: write an inverse() function for each energy function
                        mean = inverseWeightedNegLog(energy, perturbation_param_.division_weight); // convert energy to probability
                        it = tr.features.find("divProb_Var");
                        if (it == tr.features.end()) {
                            throw std::runtime_error("feature divProb_Var does not exist");
                        }
                        variance = it->second[0];
                        if (variance == 0) {
                            // do not perturb
                            LOG(logDEBUG3) << "Division: variance 0. -> do not perturb";
                            return 0.;
                        }
                        perturbed_mean = sample_with_classifier_variance(mean,variance);
                        // FIXME: the respective energy function should be used here
                        new_energy_offset = weightedNegLog(perturbed_mean, perturbation_param_.division_weight) - energy;
                        LOG(logDEBUG3) << "Division: old energy: " << energy << "; new energy offset: " << new_energy_offset;
                        return new_energy_offset;
                    case Transition:
                        if (perturbation_param_.transition_weight == 0.) {
                            return 0.;
                        }
                        // this assumes that we use the NegLog as energy function
                        // TODO: write an inverse() function for each energy function
                        mean = inverseWeightedNegLog(energy, perturbation_param_.transition_weight);
                        variance = get_transition_variance(tr,tr2);
                        if (variance == 0) {
                            // do not perturb
                            LOG(logDEBUG3) << "Transition: variance 0. -> do not perturb";
                            return 0.;
                        }
                        perturbed_mean = sample_with_classifier_variance(mean,variance);
                        // FIXME: the respective energy function should be used here
                        new_energy_offset = weightedNegLog(perturbed_mean, perturbation_param_.transition_weight) - energy;
                        LOG(logDEBUG3) << "Transition: old energy: " << energy << "; new energy offset: " << new_energy_offset;
                        return new_energy_offset;
                    default:
                        if (energyIndex >= perturbation_param_.distributionParam.size()) {
                            throw std::runtime_error("sigma is not set correctly");
                        }
                        new_energy_offset = random_normal_()*perturbation_param_.distributionParam[energyIndex];
                        LOG(logDEBUG3) << "Appearance/Disappearance: new energy offset: " << new_energy_offset;
                        return new_energy_offset;
                }
            }
        case PerturbAndMAP: //Gumbel distribution
            LOG(logDEBUG4) << "PerturbAndMAP";
            //distribution parameter: beta
            if (energyIndex >= perturbation_param_.distributionParam.size()) {
                throw std::runtime_error("sigma is not set correctly");
            }
            rand = random_uniform_();
            //throw std::runtime_error("I don't think this formula is correct; debug when needed; check whether rand>=0");
            return perturbation_param_.distributionParam[energyIndex] * log(-log(rand));
        default: //i.e. MbestCPLEX, DiverseMbest
            LOG(logDEBUG4) << "DiverseMBest/MBestCPLEX: random offset 0";
            return 0;
    }
}

size_t PerturbedInferenceModel::add_div_m_best_perturbation(marray::Marray<double>& energies,
                                                               EnergyType energy_type,
                                                               size_t factorIndex)
{
    if (perturbation_param_.distributionId == DiverseMbest){
        assert(deterministic_offset_ != NULL);
        std::vector<std::vector<size_t> >& indexlist = (*deterministic_offset_)[factorIndex];
        for (std::vector<std::vector<size_t> >::iterator index = indexlist.begin(); index != indexlist.end(); index++){
            energies(index->begin()) += perturbation_param_.distributionParam[energy_type];
        }
        factorIndex++;
    }
    return factorIndex;
}

} // namespace pgmlink

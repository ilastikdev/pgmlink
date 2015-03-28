#include "pgmlink/inferencemodel.h"

namespace pgmlink
{

InferenceModel::InferenceModel(const InferenceModel::Parameter &param):
    param_(param),
    transition_predictions_(new TransitionPredictionsMap())
{

}

InferenceModel::~InferenceModel()
{

}

void InferenceModel::use_transition_prediction_cache(InferenceModel *other)
{
    if(other == NULL)
    {
        LOG(logWARNING) << "[InferenceModel] Cannot copy cache from NULL pointer" << std::endl;
    }
    else if(!other->transition_predictions_)
    {
        LOG(logWARNING) << "[InferenceModel] Not setting empty cache" << std::endl;
    }
    else
    {
        transition_predictions_ = other->transition_predictions_;
    }
}

double InferenceModel::get_transition_prob(double distance, size_t state, double alpha) {
    double prob = exp(-distance / alpha);
    if (state == 0) {
        return 1 - prob;
    }
    return prob;
}

double InferenceModel::generateRandomOffset(EnergyType energyIndex,
                                            double energy,
                                            Traxel tr,
                                            Traxel tr2,
                                            size_t state)
{
    // unperturbed inference -> no random offset
    return 0.0;
}

bool InferenceModel::callable(boost::python::object object)
{
    return 1 == PyCallable_Check(object.ptr());
}

double InferenceModel::get_transition_probability(Traxel& tr1, Traxel& tr2, size_t state) {
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

size_t InferenceModel::add_div_m_best_perturbation(marray::Marray<double>& energies,
                                                   EnergyType energy_type,
                                                   size_t factorIndex)
{
    return factorIndex;
}



} // namespace pgmlink

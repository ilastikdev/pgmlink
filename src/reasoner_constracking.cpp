#include <algorithm>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <string.h>
#include <sstream> 
#include <memory.h>
#include <opengm/inference/lpcplex.hxx>
#include <opengm/datastructures/marray/marray.hxx>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/python.hpp>

#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/traxels.h"
#include "pgmlink/constraint_pool.hxx"

//added for view-support
#include "opengm/opengm.hxx"
#include "opengm/graphicalmodel/graphicalmodel.hxx"
#include "opengm/functions/modelviewfunction.hxx"
#include "opengm/functions/view.hxx"

//for computing inverse_sigmoid
#include <boost/math/distributions/normal.hpp>


using namespace std;

namespace pgmlink {

typedef opengm::ModelViewFunction
	<pgm::OpengmModelDeprecated::ogmGraphicalModel, marray::Marray<ValueType> >
	ViewFunctionType;

typedef opengm::LPCplex<pgm::OpengmModelDeprecated::ogmGraphicalModel,
			pgm::OpengmModelDeprecated::ogmAccumulator> cplex_optimizerHG;


ConservationTracking::~ConservationTracking() {
    	/*if (pgm_ != NULL) {
    		delete pgm_;
    		pgm_ = NULL;
    	}
       if (optimizer_ != NULL) {
          delete optimizer_;
          optimizer_ = NULL;
       }*/
}

double ConservationTracking::forbidden_cost() const {
    return forbidden_cost_;
}

void ConservationTracking::printResults(HypothesesGraph& g){
	//very verbose print of solution
	property_map<arc_active_count, HypothesesGraph::base_graph>::type& active_arcs_count =
		g.get(arc_active_count());
	property_map<node_active_count, HypothesesGraph::base_graph>::type& active_nodes_count =
		g.get(node_active_count());
	property_map<division_active_count, HypothesesGraph::base_graph>::type& active_divisions_count =
	        g.get(division_active_count());
	
	int c=0;
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
		c=0;
		for( std::vector<bool>::const_iterator i = active_arcs_count[a].begin(); i != active_arcs_count[a].end(); ++i){
			LOG(logDEBUG4) << *i;
			if (*i) {
				c++;
			}
		}
		LOG(logDEBUG4) << "total= "<<c;
	}
	
	for(std::map<HypothesesGraph::Node, size_t>::const_iterator it = app_node_map_.begin();
		it != app_node_map_.end(); ++it) {
		c=0;
		for( std::vector<long unsigned int>::const_iterator i = active_nodes_count[it->first].begin(); i != active_nodes_count[it->first].end(); ++i){
			LOG(logINFO) << *i;
			if (*i>0){
				c++;
			}
		}
		LOG(logINFO) << "total= "<<c<<std::endl;
	}
	LOG(logDEBUG4) << "division nodes "<<c<<std::endl;
		for(std::map<HypothesesGraph::Node, size_t>::const_iterator it = app_node_map_.begin();
			it != app_node_map_.end(); ++it) {
			c=0;
			for( std::vector<bool>::const_iterator i = active_divisions_count[it->first].begin(); i != active_divisions_count[it->first].end(); ++i){
				LOG(logDEBUG4) << *i<<" ";
				if (*i>0){
					c++;
				}
			}
		LOG(logDEBUG4) << "total= "<<c<<std::endl;
		}
	}

std::string ConservationTracking::get_export_filename(size_t iteration, const std::string& orig_file_name)
{
    std::stringstream export_filename;
    if(!orig_file_name.empty())
    {
        std::string orig_basename = orig_file_name;
        std::string extension = ".txt";

        // remove extension
        std::string::size_type extension_pos = orig_file_name.find_last_of(".");
        if(extension_pos != orig_file_name.npos)
        {
            extension = orig_file_name.substr(extension_pos);
            orig_basename = orig_file_name.substr(0, extension_pos);
        }

        export_filename << orig_basename << "_" << iteration << extension;
    }

    return export_filename.str();
}

void ConservationTracking::perturbedInference(HypothesesGraph& hypotheses, bool with_inference){

    reset();
    if(with_inference)
        solutions_.clear();

    pgm_ = boost::shared_ptr < pgm::OpengmModelDeprecated > (new pgm::OpengmModelDeprecated());

    HypothesesGraph *graph;

    // for formulate, add_constraints, add_finite_factors: distinguish graph & tracklet_graph
    if (with_tracklets_) {
        LOG(logINFO) << "ConservationTracking::perturbedInference: generating tracklet graph";
        tracklet2traxel_node_map_ = generateTrackletGraph2(hypotheses, tracklet_graph_);
        graph = &tracklet_graph_;
    } else {
        graph = &hypotheses;
    }

    graph->add(relative_uncertainty()).add(node_active_count());

    LOG(logINFO) << "ConservationTracking::perturbedInference: number of iterations: " << uncertainty_param_.numberOfIterations;
    LOG(logINFO) << "ConservationTracking::perturbedInference: perturb using method with Id " << uncertainty_param_.distributionId;
    LOG(logDEBUG) << "ConservationTracking::perturbedInference: formulate ";
    formulate(*graph);

    MAPGmType* model = pgm_->Model();

    LOG(logDEBUG) << "ConservationTracking::formulate: add_finite_factors";

    add_finite_factors(*graph,model,false);

    PertGmType perturbed_model = PertGmType(model->space());
    add_finite_factors(*graph, &perturbed_model, false /*perturb*/);

    LOG(logDEBUG) << "ConservationTracking::formulate: finished add_finite_factors";

    LOG(logINFO) << "number_of_transition_nodes_ = " << number_of_transition_nodes_;
    LOG(logINFO) << "number_of_appearance_nodes_ = " << number_of_appearance_nodes_;
    LOG(logINFO) << "number_of_disappearance_nodes_ = " << number_of_disappearance_nodes_;
    LOG(logINFO) << "number_of_division_nodes_ = " << number_of_division_nodes_;

    LOG(logDEBUG) <<"ConservationTracking::perturbedInference: uncertainty parameter print";

    cplex_optimizer::Parameter cplex_param;
    cplex_param.verbose_ = true;
    cplex_param.integerConstraint_ = true;
    cplex_param.epGap_ = ep_gap_;
    cplex_param.timeLimit_ = cplex_timeout_;
    uncertainty_param_.print();

    //m-best: if perturbation is set to m-best, specify number of solutions. Otherwise, we expect only one solution.
    size_t numberOfSolutions = 1;
    if (uncertainty_param_.distributionId==MbestCPLEX){
        numberOfSolutions = uncertainty_param_.numberOfIterations;
    }

    optimizer_ = new cplex_optimizer(perturbed_model,
                                     cplex_param,
                                     numberOfSolutions,
                                     get_export_filename(0, features_file_),
                                     constraints_file_,
                                     get_export_filename(0, ground_truth_file_),
                                     export_from_labeled_graph_);

    if(with_inference)
    {
        if (with_constraints_) {
            LOG(logINFO) << "add_constraints";
            add_constraints(*graph);
        }

        LOG(logINFO) << "infer MAP";
        infer();
        LOG(logINFO) << "conclude MAP";
        conclude(hypotheses);

        if(export_from_labeled_graph_ and not ground_truth_file_.empty()){
            LOG(logINFO) << "export graph labels to " << ground_truth_file_ << std::endl;
            write_labeledgraph_to_file(hypotheses);
        }

        optimizer_->set_export_file_names("","","");

        for (size_t k=1;k<numberOfSolutions;++k){
            LOG(logINFO) << "conclude "<<k+1<<"-best solution";
            optimizer_->set_export_file_names("","",get_export_filename(k, ground_truth_file_));

            opengm::InferenceTermination status = optimizer_->arg(solutions_.back(),k);
            if (status != opengm::NORMAL) {
                throw runtime_error("GraphicalModel::infer(): solution extraction terminated abnormally");
            }
            conclude(hypotheses);
        }
    }

    //store offset by factor index
    //technical assumption: the number & order of factors remains the same for the perturbed models
    //performance assumption: the number of pertubations is rather low, so iterating over
    //a list of offset for each iteration Step is not the bottle nec (this is O(n*n) time in the number of pertubations)

    size_t numberOfIterations = uncertainty_param_.numberOfIterations;
    if (uncertainty_param_.distributionId==MbestCPLEX){
        numberOfIterations = 1;
    }
    size_t num_factors = perturbed_model.numberOfFactors();
    vector<vector<vector<size_t> > >deterministic_offset(num_factors);

    //deterministic & non-deterministic perturbation
    for (size_t iterStep=1; iterStep<numberOfIterations; ++iterStep){

        num_factors = perturbed_model.numberOfFactors();
        if (uncertainty_param_.distributionId==DiverseMbest){

            for(size_t factorId=0; factorId<num_factors; ++factorId) {
                PertGmType::FactorType factor = perturbed_model[factorId];
                vector<size_t> varIndices;
                for (PertGmType::FactorType::VariablesIteratorType ind=factor.variableIndicesBegin();ind!=factor.variableIndicesEnd();++ind){
                    varIndices.push_back(solutions_[iterStep - 1][*ind]);
                }
                deterministic_offset[factorId].push_back(varIndices);
            }
        }

        LOG(logINFO) << "ConservationTracking::perturbedInference: prepare perturbation number " << iterStep;
        //initialize new model
        PertGmType perturbed_model2 = PertGmType(model->space());

        add_finite_factors(*graph, &perturbed_model2, true /*perturb*/, &deterministic_offset);
        LOG(logINFO) << "Model for perturbed inference has num factors: " << perturbed_model2.numberOfFactors();

        LOG(logINFO) << "ConservationTracking::perturbedInference construct perturbed model";
        optimizer_ = new cplex_optimizer(perturbed_model2,
                                         cplex_param,
                                         1,
                                         get_export_filename(iterStep, features_file_),
                                         "",
                                         get_export_filename(iterStep, ground_truth_file_),
                                         false);

        if(with_inference)
        {
            if (with_constraints_) {
                add_constraints(*graph);
            }
            LOG(logINFO) << "infer ";
            infer();

            LOG(logINFO) << "conclude";
            conclude(hypotheses);
        }
    }

    graph->add(relative_uncertainty());

    property_map<node_active_count, HypothesesGraph::base_graph>::type& active_nodes = graph->get(node_active_count());
    property_map<relative_uncertainty, HypothesesGraph::base_graph>::type& rel_uncertainty = graph->get(relative_uncertainty());
    for (HypothesesGraph::NodeIt n(*graph); n != lemon::INVALID; ++n) {
        double count=0;
        vector<size_t>* active_list = &active_nodes.get_value(n);
        for (vector<size_t>::iterator is_active = active_list->begin();is_active!=active_list->end();is_active++){
            if (*is_active!=0){
                ++count;
            }
        }

        rel_uncertainty.set(n,count/uncertainty_param_.numberOfIterations);
    }
}

boost::math::normal standard_gaussian_distribution(0.0, 1.0);
double sigmoid(double x){
	return cdf(standard_gaussian_distribution, x);
}

double inverse_sigmoid(double x){
	return quantile(standard_gaussian_distribution, x);
}

double ConservationTracking::sample_with_classifier_variance(double mean, double variance){
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

namespace{
    double inverseWeightedNegLog(double energy, double weight) {
        assert (weight != 0);
        return exp(-energy/weight);
    }

    double weightedNegLog(double prob, double weight) {
        if (prob <= 0.00001) {
            prob = 0.00001;
        }
        return -weight * log(prob);
    }
}
double ConservationTracking::generateRandomOffset(EnergyType energyIndex, double energy, Traxel tr, Traxel tr2, size_t state) {

    LOG(logDEBUG4) << "generateRandomOffset()";

    double rand;

    switch (uncertainty_param_.distributionId) {
		case GaussianPertubation: //normal distribution
            LOG(logDEBUG4) << "GaussianPerturbation";
            if (energyIndex >= uncertainty_param_.distributionParam.size()) {
                throw std::runtime_error("sigma is not set correctly");
            }
            return random_normal_() * uncertainty_param_.distributionParam[energyIndex];
		case ClassifierUncertainty://sample from Gaussian Distribution where variance comes from Classifier
            LOG(logDEBUG4) << "ClassifierUncertainty";
            {
                double mean, variance, perturbed_mean, new_energy_offset;
                FeatureMap::const_iterator it;

                switch (energyIndex) {
                    case Detection:
                        if (detection_weight_ == 0) {
                            return 0.;
                        }
                        // this assumes that we use the NegLog as energy function
                        // TODO: write an inverse() function for each energy function
                        mean = inverseWeightedNegLog(energy, detection_weight_); // convert energy to probability
                        it = tr.features.find("detProb_Var");
                        if (it == tr.features.end()) {
                            throw std::runtime_error("feature detProb_Var does not exist");
                        }
                        if (max_number_objects_ == 1) {
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
                        new_energy_offset = weightedNegLog(perturbed_mean, detection_weight_) - energy;
                        LOG(logDEBUG3) << "Detection: old energy: " << energy << "; new energy offset: " << new_energy_offset;
                        return new_energy_offset;
                    case Division:
                        if (division_weight_ == 0) {
                            return 0.;
                        }
                        // this assumes that we use the NegLog as energy function
                        // TODO: write an inverse() function for each energy function
                        mean = inverseWeightedNegLog(energy, division_weight_); // convert energy to probability
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
                        new_energy_offset = weightedNegLog(perturbed_mean, division_weight_) - energy;
                        LOG(logDEBUG3) << "Division: old energy: " << energy << "; new energy offset: " << new_energy_offset;
                        return new_energy_offset;
                    case Transition:
                        if (transition_weight_ == 0.) {
                            return 0.;
                        }
                        // this assumes that we use the NegLog as energy function
                        // TODO: write an inverse() function for each energy function
                        mean = inverseWeightedNegLog(energy, transition_weight_);
                        variance = get_transition_variance(tr,tr2);
                        if (variance == 0) {
                            // do not perturb
                            LOG(logDEBUG3) << "Transition: variance 0. -> do not perturb";
                            return 0.;
                        }
                        perturbed_mean = sample_with_classifier_variance(mean,variance);
                        // FIXME: the respective energy function should be used here
                        new_energy_offset = weightedNegLog(perturbed_mean, transition_weight_) - energy;
                        LOG(logDEBUG3) << "Transition: old energy: " << energy << "; new energy offset: " << new_energy_offset;
                        return new_energy_offset;
                    default:
                        if (energyIndex >= uncertainty_param_.distributionParam.size()) {
                            throw std::runtime_error("sigma is not set correctly");
                        }
                        new_energy_offset = random_normal_()*uncertainty_param_.distributionParam[energyIndex];
                        LOG(logDEBUG3) << "Appearance/Disappearance: new energy offset: " << new_energy_offset;
                        return new_energy_offset;
                }
			}
		case PerturbAndMAP: //Gumbel distribution
            LOG(logDEBUG4) << "PerturbAndMAP";
			//distribution parameter: beta
            if (energyIndex >= uncertainty_param_.distributionParam.size()) {
                throw std::runtime_error("sigma is not set correctly");
            }
            rand = random_uniform_();
            throw std::runtime_error("I don't think this formula is correct; debug when needed; check whether rand>=0");
            return uncertainty_param_.distributionParam[energyIndex] * log(-log(rand));
		default: //i.e. MbestCPLEX, DiverseMbest
            LOG(logDEBUG4) << "DiverseMBest/MBestCPLEX: random offset 0";
			return 0;
	}
}

void ConservationTracking::formulate(const HypothesesGraph& hypotheses) {
    LOG(logDEBUG) << "ConservationTracking::formulate: entered";

    LOG(logDEBUG) << "ConservationTracking::formulate: add_transition_nodes";
    add_transition_nodes(hypotheses);
    LOG(logDEBUG) << "ConservationTracking::formulate: add_appearance_nodes";
    add_appearance_nodes(hypotheses);
    LOG(logDEBUG) << "ConservationTracking::formulate: add_disappearance_nodes";
    add_disappearance_nodes(hypotheses);

    LOG(logDEBUG) << "ConservationTracking::formulate: add_division_nodes";
    if (with_divisions_) {
        add_division_nodes(hypotheses);
    }
}

void ConservationTracking::infer() {
	if (!with_constraints_) {
		//opengm::hdf5::save(optimizer_->graphicalModel(), "./conservationTracking.h5", "conservationTracking");
		throw std::runtime_error("GraphicalModel::infer(): inference with soft constraints is not implemented yet. The conservation tracking factor graph has been saved to file");
	}
    opengm::InferenceTermination status = optimizer_->infer();
    if (status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): optimizer terminated abnormally");
    }

    solutions_.push_back(IlpSolution());
    opengm::InferenceTermination statusExtract = optimizer_->arg(solutions_.back());
    if (statusExtract != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): solution extraction terminated abnormally");
    }
}

void ConservationTracking::conclude( HypothesesGraph& g) {
    if(export_from_labeled_graph_ and not  ground_truth_file_.empty()){
        clpex_variable_id_map_ = optimizer_->get_clpex_variable_id_map();
        clpex_factor_id_map_ = optimizer_->get_clpex_factor_id_map();
    }

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

    if (!with_tracklets_) {
        tracklet_graph_.add(tracklet_intern_arc_ids()).add(traxel_arc_id());
    }
    property_map<tracklet_intern_arc_ids, HypothesesGraph::base_graph>::type& tracklet_arc_id_map =
            tracklet_graph_.get(tracklet_intern_arc_ids());
    property_map<traxel_arc_id, HypothesesGraph::base_graph>::type& traxel_arc_id_map =
            tracklet_graph_.get(traxel_arc_id());

	int iterStep = active_nodes_count[app_node_map_.begin()->first].size();
	bool isMAP = (iterStep==0);

	if (isMAP){
		//initialize vectors for storing optimizer results
		for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
		active_arcs_count.set(a,std::vector<bool>());
		}
		for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
			active_nodes_count.set(n,std::vector<long unsigned int>());
			active_divisions_count.set(n, std::vector<bool>());
		}
	}

	//initialize node counts by 0
	for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
			active_nodes_count.get_value(n).push_back(0);
			active_divisions_count.get_value(n).push_back(0);
			}


	//initialize arc counts by 0
	for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
		active_arcs.set(a, false);
		active_arcs_count.get_value(a).push_back(0);
	}
    // write state after inference into 'active'-property maps
    // the node is also active if its appearance node is active
    for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = app_node_map_.begin();
            it != app_node_map_.end(); ++it) {
        if (with_tracklets_) {
            // set state of tracklet nodes
            std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map_[it->first];

            for (std::vector<HypothesesGraph::Node>::const_iterator tr_n_it = traxel_nodes.begin();
                    tr_n_it != traxel_nodes.end(); ++tr_n_it) {
                HypothesesGraph::Node n = *tr_n_it;
                active_nodes.set(n, solutions_.back()[it->second]);
                active_nodes_count.get_value(n)[iterStep]=solutions_.back()[it->second];
				//TODO: active_nodes_vector
            }

            // set state of tracklet internal arcs
            std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
            for (std::vector<int>::const_iterator arc_id_it = arc_ids.begin();
                    arc_id_it != arc_ids.end(); ++arc_id_it) {
                HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
                assert(active_arcs[a] == false);
                if (solutions_.back()[it->second] > 0) {

					active_arcs.set(a, true);
					active_arcs_count.get_value(a)[iterStep]=true;

                    assert(active_nodes[g.source(a)] == solutions_.back()[it->second]
                            && "tracklet internal arcs must have the same flow as their connected nodes");
                    assert(active_nodes[g.target(a)] == solutions_.back()[it->second]
                            && "tracklet internal arcs must have the same flow as their connected nodes");
                }
            }
        } else {

            active_nodes.set(it->first, solutions_.back()[it->second]);
            active_nodes_count.get_value(it->first)[iterStep]=solutions_.back()[it->second];
        }
    }
    // the node is also active if its disappearance node is active
    for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = dis_node_map_.begin();
            it != dis_node_map_.end(); ++it) {

        if (solutions_.back()[it->second] > 0) {
            if (with_tracklets_) {
                // set state of tracklet nodes
                std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map_[it->first];
                for (std::vector<HypothesesGraph::Node>::const_iterator tr_n_it =
                        traxel_nodes.begin(); tr_n_it != traxel_nodes.end(); ++tr_n_it) {
                    HypothesesGraph::Node n = *tr_n_it;

                    if (active_nodes[n] == 0) {

                        active_nodes.set(n, solutions_.back()[it->second]);
                        active_nodes_count.get_value(n)[iterStep]=solutions_.back()[it->second];

                    } else {
                        assert(active_nodes[n] == solutions_.back()[it->second]);
                    }
                }
                // set state of tracklet internal arcs
                std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
                for (std::vector<int>::const_iterator arc_id_it = arc_ids.begin();
                        arc_id_it != arc_ids.end(); ++arc_id_it) {
                    HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
                    if (solutions_.back()[it->second] > 0) {

						active_arcs.set(a, true);
						active_arcs_count.get_value(a)[iterStep]=true;

                        assert(active_nodes[g.source(a)] == solutions_.back()[it->second]
                                && "tracklet internal arcs must have the same flow as their connected nodes");
                        assert(active_nodes[g.target(a)] == solutions_.back()[it->second]
                                && "tracklet internal arcs must have the same flow as their connected nodes");
                    }
                }
            } else {

				if (active_nodes[it->first] == 0) {
                    active_nodes.set(it->first, solutions_.back()[it->second]);
                    active_nodes_count.get_value(it->first)[iterStep]=solutions_.back()[it->second];

                } else{
                    assert(active_nodes[it->first] == solutions_.back()[it->second]);
                }
            }
        }
    }

    for (std::map<HypothesesGraph::Arc, size_t>::const_iterator it = arc_map_.begin();
            it != arc_map_.end(); ++it) {
        if (solutions_.back()[it->second] >= 1) {
            if (with_tracklets_) {
				active_arcs.set(g.arcFromId((traxel_arc_id_map[it->first])), true);
				active_arcs_count.get_value(g.arcFromId((traxel_arc_id_map[it->first])))[iterStep]=true;
            } else {
				active_arcs.set(it->first, true);
				active_arcs_count.get_value(it->first)[iterStep]=true;
            }
        }
    }
    // write division node map
    if (with_divisions_) {
    	for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = div_node_map_.begin();
    	                it != div_node_map_.end(); ++it) {
    	            division_nodes.set(it->first, false);
    	        }
        for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = div_node_map_.begin();
                it != div_node_map_.end(); ++it) {

            if (solutions_.back()[it->second] >= 1) {
                if (with_tracklets_) {
                    // set division property for the last node in the tracklet
					HypothesesGraph::Node n = tracklet2traxel_node_map_[it->first].back();
                    division_nodes.set(n, true);
                    active_divisions_count.get_value(n)[iterStep]=true;
                } else {
                    division_nodes.set(it->first, true);
                    active_divisions_count.get_value(it->first)[iterStep]=true;
                }
            }
        }
	}
}

const std::map<HypothesesGraph::Arc, size_t>& ConservationTracking::get_arc_map() const {
    return arc_map_;
}

const std::vector<ConservationTracking::IlpSolution> &ConservationTracking::get_ilp_solutions() const
{
    return solutions_;
}

void ConservationTracking::set_ilp_solutions(const std::vector<ConservationTracking::IlpSolution>& solutions)
{
    if(solutions.size() == 0)
    {
        LOG(logWARNING) << "No solutions given to set";
        return;
    }

    solutions_.clear();
    for(std::vector<ConservationTracking::IlpSolution>::const_iterator sol_it = solutions.begin();
        sol_it != solutions.end();
        ++sol_it)
    {
        solutions_.push_back(ConservationTracking::IlpSolution(*sol_it));
    }
}

void ConservationTracking::reset() {
    if (optimizer_ != NULL) {
        delete optimizer_;
        optimizer_ = NULL;
    }
    arc_map_.clear();
    div_node_map_.clear();
    app_node_map_.clear();
    dis_node_map_.clear();
//    solutions_.clear();
}

void ConservationTracking::add_appearance_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        pgm_->Model()->addVariable(max_number_objects_ + 1);
        app_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;

        HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
        nodes_per_timestep_[timestep_map[n]].push_back(pgm_->Model()->numberOfVariables() - 1);

        assert(pgm_->Model()->numberOfLabels(app_node_map_[n]) == max_number_objects_ + 1);
        ++count;
    }
    number_of_appearance_nodes_ = count;
}

void ConservationTracking::add_disappearance_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        pgm_->Model()->addVariable(max_number_objects_ + 1);
        dis_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;

        HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
        nodes_per_timestep_[timestep_map[n]].push_back(pgm_->Model()->numberOfVariables() - 1);

        assert(pgm_->Model()->numberOfLabels(dis_node_map_[n]) == max_number_objects_ + 1);
        ++count;
    }
    number_of_disappearance_nodes_ = count;
}

void ConservationTracking::add_transition_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        pgm_->Model()->addVariable(max_number_objects_ + 1);
        arc_map_[a] = pgm_->Model()->numberOfVariables() - 1;
        // store these nodes by the timestep of the base-appearance node,
        // as well as in the timestep of the disappearance node they enter
        HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
        HypothesesGraph::OutArcIt out_arc(g, a);
        HypothesesGraph::Node n = g.source(out_arc);
        nodes_per_timestep_[timestep_map[n]].push_back(pgm_->Model()->numberOfVariables() - 1);
        n = g.target(out_arc);
        nodes_per_timestep_[timestep_map[n]].push_back(pgm_->Model()->numberOfVariables() - 1);

        assert(pgm_->Model()->numberOfLabels(arc_map_[a]) == max_number_objects_ + 1);
        ++count;
    }
    number_of_transition_nodes_ = count;
}

void ConservationTracking::add_division_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        size_t number_of_outarcs = 0;
        for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
            ++number_of_outarcs;
        }
        if (number_of_outarcs > 1) {
            pgm_->Model()->addVariable(2);
            div_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
            HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
            nodes_per_timestep_[timestep_map[n]].push_back(pgm_->Model()->numberOfVariables() - 1);

            assert(pgm_->Model()->numberOfLabels(div_node_map_[n]) == 2);
            ++count;
        }
    }
    number_of_division_nodes_ = count;
}

namespace {
double get_transition_prob(double distance, size_t state, double alpha) {
    double prob = exp(-distance / alpha);
    if (state == 0) {
        return 1 - prob;
    }
    return prob;
}
}

boost::python::dict convertFeatureMapToPyDict(FeatureMap map){
	boost::python::dict dictionary;
	for (FeatureMap::iterator iter = map.begin(); iter != map.end(); ++iter) {
			dictionary[iter->first] = iter->second;
		}
	return dictionary;
}


double ConservationTracking::get_transition_variance(Traxel tr1, Traxel tr2) {
    double var;

    if (transition_classifier_.ptr()==boost::python::object().ptr()){
        var = uncertainty_param_.distributionParam[Transition];
        LOG(logDEBUG4) << "using constant transition variance " << var;
        if (var < 0) {
            throw std::runtime_error("the transition variance must be positive");
        }
    } else {        
        TransitionPredictionsMap::const_iterator it = transition_predictions_.find(std::make_pair(tr1, tr2));
        if ( it == transition_predictions_.end() ) {
            throw std::runtime_error("cannot find prob/var. get_transition_probability must be called first");
        }
        var = it->second.second;
        LOG(logDEBUG4) << "using GPC transition variance " << var;
    }

    return var;
}

double ConservationTracking::get_transition_probability(Traxel tr1, Traxel tr2, size_t state) {

    LOG(logDEBUG4) << "get_transition_probability()";

    double prob;

    //read the FeatureMaps from Traxels
    if (transition_classifier_.ptr()==boost::python::object().ptr()) {
        double distance = 0;
        if (with_optical_correction_) {
            distance = tr1.distance_to_corr(tr2);
        } else {
            distance = tr1.distance_to(tr2);
        }
        prob = get_transition_prob(distance, state, transition_parameter_);
        LOG(logDEBUG4) << "get_transition_probability(): using deterministic function: " << tr1
                       << " " << tr2 << " [" << state << "] = " << prob << "; distance = " << distance;
        assert(prob >= 0 && prob <= 1);
        return prob;
    }

    TransitionPredictionsMap::const_iterator it = transition_predictions_.find(std::make_pair(tr1, tr2));
    if ( it == transition_predictions_.end() ) {
        // predict and store
        double var;
        try {            
            boost::python::object prediction = transition_classifier_.attr("predict")(tr1,tr2);

            boost::python::object probs_python = prediction.attr("__getitem__")(0);
            // we are only interested in the probability of the second class, since it is a binary classifier
            prob = boost::python::extract<double>(probs_python.attr("__getitem__")(1));
            var = boost::python::extract<double>(prediction.attr("__getitem__")(1));

            LOG(logDEBUG4) << "python::transition_classifier, prob =" << prob << "; var =" << var;
        } catch (...) {
            throw std::runtime_error("cannot call the transition classifier from python");
        }
        transition_predictions_[std::make_pair(tr1, tr2)] = std::make_pair(prob, var);
    } else {
        prob = it->second.first;
    }

	if (state == 0) {
        prob = 1-prob;
    }
    LOG(logDEBUG4) << "get_transition_probability(): using Gaussian process classifier: p[" << state << "] = " << prob ;
    return prob;
}

template <typename ModelType>
void ConservationTracking::add_finite_factors(const HypothesesGraph& g, ModelType* model, bool perturb /*=false*/, vector< vector<vector<size_t> > >* detoff /*= NULL*/) {
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

	size_t factorIndex = 0;
    if(model->numberOfFactors()!=factorIndex)
        LOG(logERROR) << "Model has a too large number of factors: " << model->numberOfFactors() << ">" << factorIndex;
	assert(model->numberOfFactors()==factorIndex);

    LOG(logDEBUG) << "ConservationTracking::add_finite_factors: entered";
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
    		g.get(node_tracklet());

//    bool perturb_transitions_locally=(perturb && param_.distributionId==ClassifierUncertainty);
//      //if transitions ought to be perturbed, generate offset for RegionCenters in order to perturb distances->probabilities->energies
//    map<Traxel,vector<double> > offset;
//    if (perturb_transitions_locally){
//    	Traxel tr;
//    	for (HypothesesGraph::NodeIt it(g); it != lemon::INVALID; ++it) {

//    		if (with_tracklets_) {
//    			std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map_[it];
//    			for (std::vector<HypothesesGraph::Node>::const_iterator tr_n_it = traxel_nodes.begin();
//    					tr_n_it != traxel_nodes.end(); ++tr_n_it) {
//    				HypothesesGraph::Node n = *tr_n_it;
//    				tr = traxel_map_[n];
//    				for (size_t dim=0;dim<tr.features["com"].size();dim++){
//                        offset[tr].push_back(generateRandomOffset(Transition));
//    				}
//    			}
//    		} else {
//    			tr = traxel_map_[it];
//    			for (size_t dim=0;dim<tr.features["com"].size();dim++){
//    				offset[tr].push_back(generateRandomOffset(Transition));
//    			}
//    		}
//    	}
//    }


    ////
    //// add detection factors
    ////
    LOG(logDEBUG) << "ConservationTracking::add_finite_factors: add detection factors";
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        size_t num_vars = 0;
        vector<size_t> vi;
        vector<double> cost;

        int node_begin_time = -1;
        int node_end_time = -1;
        if (with_tracklets_) {
            node_begin_time = tracklet_map_[n].front().Timestep;
            node_end_time = tracklet_map_[n].back().Timestep;
        } else {
            node_begin_time = traxel_map_[n].Timestep;
            node_end_time = traxel_map_[n].Timestep;
        }


        double energy,e;
        if (app_node_map_.count(n) > 0) {
            vi.push_back(app_node_map_[n]);
            if (node_begin_time <= g.earliest_timestep()) {  // "<" holds if there are only tracklets in the first frame
                // pay no appearance costs in the first timestep
                cost.push_back(0.);
            } else {
            	if (with_tracklets_) {
            		energy = appearance_cost_(tracklet_map_[n].front());
            	} else {
            		energy = appearance_cost_(traxel_map_[n]);
            	}
            	if (perturb){
                    energy+= generateRandomOffset(Appearance);
            	}
            	cost.push_back(energy);
            }
            ++num_vars;
        }
		
		if (dis_node_map_.count(n) > 0) {
            vi.push_back(dis_node_map_[n]);
            if (node_end_time < g.latest_timestep()) { // "<" holds if there are only tracklets in the last frame
            	if (with_tracklets_) {
            		energy = disappearance_cost_(tracklet_map_[n].back());
            	} else {
            		energy = disappearance_cost_(traxel_map_[n]);
            	}
            	if (perturb){
                    energy += generateRandomOffset(Disappearance);
            	}
            	cost.push_back(energy);
            } else {
            	cost.push_back(0);
            }
            ++num_vars;
        }

        // convert vector to array
        vector<size_t> coords(num_vars, 0); // number of variables
        // ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
        //pgm::OpengmExplicitFactor<double> table(vi.begin(), vi.end(), forbidden_cost_, (max_number_objects_ + 1));
        vector<size_t> shape(num_vars,(max_number_objects_+1));
        marray::Marray<double> energies(shape.begin(),shape.end(),forbidden_cost_);

        for (size_t state = 0; state <= max_number_objects_; ++state) {
        	if (with_tracklets_) {
        		energy=0;
        		// add all detection factors of the internal nodes
        		for (std::vector<Traxel>::const_iterator trax_it = tracklet_map_[n].begin();
        				trax_it != tracklet_map_[n].end(); ++trax_it) {                    
        			e = detection_(*trax_it, state);
        			energy += e;                    
        			if (perturb){                        
                        energy+= generateRandomOffset(Detection, e, *trax_it, state);
        			}
        		}
        		// add all transition factors of the internal arcs                
                Traxel tr_prev;
                bool first = true;
                for (std::vector<Traxel>::const_iterator trax_it = tracklet_map_[n].begin();
                        trax_it != tracklet_map_[n].end(); ++trax_it) {
                    LOG(logDEBUG4) << "internal arcs traxel " << *trax_it;
                    Traxel tr = *trax_it;
                    if (!first) {
                        e = transition_( get_transition_probability(tr_prev, tr, state) );
                        energy += e;
                        if (perturb){
                            energy += generateRandomOffset(Transition, e, tr_prev, tr);
                        }
                    } else {
                        first = false;
                    }
                    tr_prev = tr;
                }

        	} else {                
        		e = detection_(traxel_map_[n], state);
        		energy=e;                
    			if (perturb){
                    energy+= generateRandomOffset(Detection, e, traxel_map_[n], state);
    			}
        	}            
            LOG(logDEBUG2) << "ConservationTracking::add_finite_factors: detection[" << state
            	 << "] = " << energy;
            for (size_t var_idx = 0; var_idx < num_vars; ++var_idx) {
                coords[var_idx] = state;   
                // if only one of the variables is > 0, then it is an appearance in this time frame
                // or a disappearance in the next timeframe. Hence, add the cost of appearance/disappearance
                // to the detection cost
                energies(coords.begin()) = energy + state * cost[var_idx];
                coords[var_idx] = 0;
                LOG(logDEBUG4) << "ConservationTracking::add_finite_factors: var_idx "
                                << var_idx << " = " << energy;
            }
            // also this energy if both variables have the same state
            if (num_vars == 2) {
                coords[0] = state;
                coords[1] = state;
                // only pay detection energy if both variables are on
                //table.set_value(coords, energy);
                energies(coords.begin())=energy;
                coords[0] = 0;
                coords[1] = 0;

                LOG(logDEBUG4) << "ConservationTracking::add_finite_factors: var_idxs 0 and var_idx 1 = "
                                << energy;
            }
        }

        LOG(logDEBUG3) << "ConservationTracking::add_finite_factors: adding table to pgm";
        //functor add detection table
        if (perturb && uncertainty_param_.distributionId==DiverseMbest){
        	vector<vector<size_t> >* indexlist = &detoff->operator[](factorIndex);
        	for (vector<vector<size_t> >::iterator index=indexlist->begin();index !=indexlist->end();index++){
                energies(index->begin()) += uncertainty_param_.distributionParam[Detection];
        	}
        	factorIndex++;
        }
        typename ModelType::FunctionIdentifier funcId = model->addFunction(energies);

        sort(vi.begin(),vi.end());
        model->addFactor(funcId,vi.begin(),vi.end());
        if (not perturb)
            detection_f_node_map_[n] = model->numberOfFactors() -1;
        //table.add_to(model);
    }

    ////
    //// add transition factors
    ////
    LOG(logDEBUG) << "ConservationTracking::add_finite_factors: add transition factors";

    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        size_t vi[] = { arc_map_[a] };
        vector<size_t> coords(1, 0); // number of variables
        // ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
        // pgm::OpengmExplicitFactor<double> table(vi, vi + 1, forbidden_cost_, (max_number_objects_ + 1));

        vector<size_t> shape(1,(max_number_objects_+1));
        marray::Marray<double> energies(shape.begin(),shape.end(),forbidden_cost_);

        for (size_t state = 0; state <= max_number_objects_; ++state) {
            Traxel tr1, tr2;
            if (with_tracklets_) {
                tr1 = tracklet_map_[g.source(a)].back();
                tr2 = tracklet_map_[g.target(a)].front();
            } else {
                tr1 = traxel_map_[g.source(a)];
                tr2 = traxel_map_[g.target(a)];
            }

            double energy = transition_(get_transition_probability(tr1, tr2, state));

            if (perturb) {
                energy += generateRandomOffset(Transition, energy, tr1, tr2);
        	}

        	LOG(logDEBUG2) << "ConservationTracking::add_finite_factors: transition[" << state
        			<< "] = " << energy;
        	coords[0] = state;
        	//table.set_value(coords, energy);
        	energies(coords.begin())=energy;
        	coords[0] = 0;
        }
        if (perturb && uncertainty_param_.distributionId==DiverseMbest){
        	vector<vector<size_t> >* indexlist = &detoff->operator[](factorIndex);
        	for (vector<vector<size_t> >::iterator index=indexlist->begin();index !=indexlist->end();index++){
                energies(index->begin())+=uncertainty_param_.distributionParam[Transition];
        	}
        	factorIndex++;
        }
        typename ModelType::FunctionIdentifier funcId = model->addFunction(energies);
        model->addFactor(funcId,vi,vi+1);
    }
       
    ////
    //// add division factors
    ////
    if (with_divisions_) {
        LOG(logDEBUG) << "ConservationTracking::add_finite_factors: add division factors";
        for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
            if (div_node_map_.count(n) == 0) {
                continue;
            }
            size_t vi[] = { div_node_map_[n] };
            vector<size_t> coords(1, 0); // number of variables
            // ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
            //pgm::OpengmExplicitFactor<double> table(vi, vi + 1, forbidden_cost_, 2);
            vector<size_t> shape(1,2);
            marray::Marray<double> energies(shape.begin(),shape.end(),forbidden_cost_);

            for (size_t state = 0; state <= 1; ++state) {
            	double energy;
            	Traxel tr;
            	if (with_tracklets_) {
            		tr = tracklet_map_[n].back();
            	} else {
            		tr = traxel_map_[n];
            	}
        		energy = division_(tr, state);
            	if (perturb){
            		energy+= generateRandomOffset(Division, energy,  tr);
            	}
                LOG(logDEBUG2) << "ConservationTracking::add_finite_factors: division[" << state
                        << "] = " << energy;
                coords[0] = state;
                //table.set_value(coords, energy);
                energies(coords.begin())=energy;
                coords[0] = 0;
            }
            //table.add_to(model);
            if (perturb && uncertainty_param_.distributionId==DiverseMbest){
            	vector<vector<size_t> >* indexlist = &detoff->operator[](factorIndex);
            	for (vector<vector<size_t> >::iterator index=indexlist->begin();index !=indexlist->end();index++){
                    energies(index->begin())+=uncertainty_param_.distributionParam[Division];
            	}
            	factorIndex++;
            }
            typename ModelType::FunctionIdentifier funcId = model->addFunction(energies);
            model->addFactor(funcId,vi,vi+1);
        }
    }
}

void ConservationTracking::write_labeledgraph_to_file(const HypothesesGraph& g){

    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    property_map<traxel_arc_id, HypothesesGraph::base_graph>::type& traxel_arc_id_map = tracklet_graph_.get(traxel_arc_id());

    //write this map to label file
    map<int,label_type> cplexid_label_map;
    map<int,stringstream> cplexid_label_info_map;
    map<size_t,size_t> cplexid_weight_class_map;

    // fill labels of Variables
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        //appearance 
        for (size_t state = 0; state < pgm_->Model()->numberOfLabels(app_node_map_[n]); ++state) {
            int id = clpex_variable_id_map_[make_pair(app_node_map_[n],state)];
            cplexid_label_map[id] = ((g.get(appearance_label())[n]==state)?1:0);
            LOG(logDEBUG4) <<"app\t"<< cplexid_label_map[id] << "  " << id << "  " <<  app_node_map_[n] << "  " << state;
            cplexid_label_info_map[id] <<  "# appearance id:" << id << " state:" << g.get(appearance_label())[n] << "/" << state << "  node:" << app_node_map_[n]  << "traxel id:" <<  traxel_map[n].Id << "  ts:" << traxel_map[n].Timestep ;
            cplexid_weight_class_map[id] = 0;     
        }    
        //disappearance
        for (size_t state = 0; state < pgm_->Model()->numberOfLabels(dis_node_map_[n]); ++state) { 
            int id = clpex_variable_id_map_[make_pair(dis_node_map_[n],state)];
            cplexid_label_map[id] = ((g.get(disappearance_label())[n]==state)?1:0);
            LOG(logDEBUG4) <<"dis\t"<< cplexid_label_map[id] << "  " << id << "  " <<  dis_node_map_[n] << "  " << state;
            cplexid_label_info_map[id] <<  "# disappearance id:" << id << " state:" << g.get(disappearance_label())[n] << "/" << state << "  node:" << dis_node_map_[n]  << "traxel id:" <<  traxel_map[n].Id << "  ts:" << traxel_map[n].Timestep ;
            cplexid_weight_class_map[id] = 1;
        }    
        //division
        if(with_divisions_ and div_node_map_.count(n) != 0){
                
            for (size_t state = 0; state < pgm_->Model()->numberOfLabels(div_node_map_[n]); ++state) {
                int id = clpex_variable_id_map_[make_pair(div_node_map_[n],state)];
                cplexid_label_map[id] = ((g.get(division_label())[n]==state)?1:0);
                LOG(logDEBUG4) <<"div\t"<< cplexid_label_map[id] << "  " << id << "  " <<  div_node_map_[n] << "  " << state <<"   "<< number_of_division_nodes_;
                cplexid_label_info_map[id] <<  "# division id:" << id << " state:" << g.get(division_label())[n] << "/" << state << "  node:" << div_node_map_[n]  << "traxel id:" <<  traxel_map[n].Id << "  ts:" << traxel_map[n].Timestep ;
                cplexid_weight_class_map[id] = 2;
            }    
        }
    }
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        //move
        for (size_t state = 0; state < pgm_->Model()->numberOfLabels(arc_map_[a]); ++state) {
            int id = clpex_variable_id_map_[make_pair(arc_map_[a],state)];
            cplexid_label_map[id] = ((g.get(arc_label())[a]==state)?1:0);
            LOG(logDEBUG4) <<"arc\t"<< cplexid_label_map[id] << "  " <<id << "  " <<  arc_map_[a] << "  " << state;
            cplexid_label_info_map[id] <<  "# move id:" << id << " state:" << g.get(arc_label())[a] << "/" << state << "  node:" << arc_map_[a]  << "traxel id:" <<  traxel_map[g.source(a)].Id << "-->" << traxel_map[g.target(a)].Id << "  ts:" << traxel_map[g.target(a)].Timestep ;
            cplexid_weight_class_map[id] = 3;
            //traxel_arc_id_map[a].Timestep
        }
    }

    // fill labels of Factors (only second order factors need to be exported (others accounted for in variable states))
    
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        //detection factor detection_node_map_
        for (size_t s1 = 0; s1 < pgm_->Model()->numberOfLabels(detection_f_node_map_[n]); ++s1) {
            for (size_t s2 = 0; s2 < pgm_->Model()->numberOfLabels(detection_f_node_map_[n]); ++s2) {
                int id = clpex_factor_id_map_[make_pair(detection_f_node_map_[n],make_pair(s1,s2))];
                cplexid_label_map[id] = ((g.get(appearance_label())[n]==s1 and g.get(disappearance_label())[n]==s2)?1:0);
                LOG(logDEBUG4) <<"detection\t"<< cplexid_label_map[id] <<"  "<<id<< "  " <<  detection_f_node_map_[n] << "  " << s1 <<"  " << s2 << endl;
                cplexid_weight_class_map[id] = 4;
                cplexid_label_info_map[id] <<  "# factor id:" << id << " state:" << g.get(appearance_label())[n] <<","<<g.get(disappearance_label())[n]<< "/" << s1 << "," << s2 << "  node:" << app_node_map_[n]  << "traxel id:" <<  traxel_map[n].Id << "  ts:" << traxel_map[n].Timestep ;
            }
        }
    }

    

    //write labels to file
    std::ofstream ground_truth_file;
        
    ground_truth_file.open (ground_truth_file_,std::ios::app);
        
    for(std::map<int,size_t>::iterator iterator = cplexid_label_map.begin(); iterator != cplexid_label_map.end(); iterator++) {
        ground_truth_file << iterator->second <<"\t\t#c"<<cplexid_weight_class_map[iterator->first]<<"\t\tcplexid:" << iterator->first  << cplexid_label_info_map[iterator->first].str() <<endl;
    }
    ground_truth_file.close();

}

size_t ConservationTracking::cplex_id(size_t opengm_id, size_t state) {
    return optimizer_->lpNodeVi(opengm_id, state);
}

//set up optimizer from constraints by reading from formulated gm
void ConservationTracking::add_constraints(const HypothesesGraph& g) {
    LOG(logDEBUG) << "ConservationTracking::add_constraints: entered";

    pgm::ConstraintPool constraint_pool(forbidden_cost_, with_divisions_, with_appearance_, with_disappearance_, with_misdetections_allowed_);

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        ////
        //// outgoing transitions
        ////
        {
            std::vector<size_t> transition_nodes;
            for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
                transition_nodes.push_back(arc_map_[a]);
            }

            int division_node = -1;
            if(div_node_map_.count(n) > 0)
            {
                division_node = div_node_map_[n];
            }
            size_t appearance_node = app_node_map_[n];

            constraint_pool.add_constraint(pgm::ConstraintPool::OutgoingConstraint(appearance_node, division_node, transition_nodes));
        }


        ////
        //// incoming transitions
        ////
        {
            std::vector<size_t> transition_nodes;
            for (HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a) {
                transition_nodes.push_back(arc_map_[a]);
            }
            size_t disappearance_node = dis_node_map_[n];

            constraint_pool.add_constraint(pgm::ConstraintPool::IncomingConstraint(transition_nodes, disappearance_node));
        }

        ////
        //// disappearance/appearance coupling
        ////
        if (app_node_map_.count(n) > 0 && dis_node_map_.count(n) > 0)
        {
            constraint_pool.add_constraint(pgm::ConstraintPool::DetectionConstraint((size_t)dis_node_map_[n], (size_t)app_node_map_[n]));
        }
    }

    constraint_pool.force_softconstraint(!with_constraints_);


//    // For testing purposes:
//    boost::posix_time::time_facet *facet = new boost::posix_time::time_facet("%Y%m%d-%H%M%S");
//    std::stringstream filename_model;
//    filename_model.imbue(std::locale(std::cout.getloc(), facet));
//    filename_model << "./constracking_model-" << boost::posix_time::second_clock::local_time() << ".h5";

//    opengm::hdf5::save(*pgm_->Model(), filename_model.str(), "model");
    
//    std::stringstream filename_constraints;
//    filename_constraints.imbue(std::locale(std::cout.getloc(), facet));
//    filename_constraints << "./constracking_constraints-" << boost::posix_time::second_clock::local_time() << ".cp";

//    std::ofstream out_stream(filename_constraints.str().c_str());
//    boost::archive::text_oarchive oa(out_stream);
//    oa & constraint_pool;
//    oa & nodes_per_timestep_;
//    out_stream.close();

    constraint_pool.add_constraints_to_problem(*pgm_->Model(), *optimizer_);
}

} /* namespace pgmlink */

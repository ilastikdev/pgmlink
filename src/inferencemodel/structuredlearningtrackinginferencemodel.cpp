//#include "pgmlink/inferencemodel/constrackingexplicitinferencemodel.h"
#include <opengm/functions/learnable/lsum_of_experts.hxx>
#include "pgmlink/inferencemodel/structuredlearningtrackinginferencemodel.h"
#include <boost/python.hpp>

//#include "opengm/opengm.hxx"

namespace pgmlink
{

//void StructuredLearningTrackingInferenceModel::setWeight(size_t index, double val){
//    weights_[index] = val;
//    std::cout << " ===================================================" << weights_[index] << std::endl;
//}
//double StructuredLearningTrackingInferenceModel::weight(size_t index){
//    std::cout << " ===================================================" << weights_[index] << std::endl;
//    return weights_[index];
//}


size_t StructuredLearningTrackingInferenceModel::add_detection_factors(const HypothesesGraph& g, size_t factorIndex)
{
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_detection_factors" << weights_.getWeight((size_t) 0) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_detection_factors" << weights_.getWeight((size_t) 1) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_detection_factors" << weights_.getWeight((size_t) 2) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_detection_factors" << weights_.getWeight((size_t) 3) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_detection_factors" << weights_.getWeight((size_t) 4) << std::endl;
    ////
    //// add detection factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors" << std::endl;
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());
    //std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors" << std::endl;

    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_finite_factors: add detection factors";
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        size_t num_vars = 0;
        std::vector<size_t> vi;
        std::vector<double> cost;

        int node_begin_time = -1;
        int node_end_time = -1;
        if (param_.with_tracklets)
        {
            node_begin_time = tracklet_map_[n].front().Timestep;
            node_end_time = tracklet_map_[n].back().Timestep;
        }
        else
        {
            node_begin_time = traxel_map_[n].Timestep;
            node_end_time = traxel_map_[n].Timestep;
        }


        double energy, e, f, w;
        if (app_node_map_.count(n) > 0)
        {
            vi.push_back(app_node_map_[n]);
            // "<" holds if there are only tracklets in the first frame
            if (node_begin_time <= g.earliest_timestep())
            {
                // pay no appearance costs in the first timestep
                cost.push_back(0.);
            }
            else
            {
                if (param_.with_tracklets)
                {
                    energy = param_.appearance_cost(tracklet_map_[n].front());
                }
                else
                {
                    energy = param_.appearance_cost(traxel_map_[n]);
                }

                energy += generateRandomOffset(Appearance);
                cost.push_back(energy);
            }
            ++num_vars;
        }

        if (dis_node_map_.count(n) > 0)
        {
            vi.push_back(dis_node_map_[n]);
            // "<" holds if there are only tracklets in the last frame
            if (node_end_time < g.latest_timestep())
            {
                if (param_.with_tracklets)
                {
                    energy = param_.disappearance_cost(tracklet_map_[n].back());
                }
                else
                {
                    energy = param_.disappearance_cost(traxel_map_[n]);
                }

                energy += generateRandomOffset(Disappearance);
                cost.push_back(energy);
            }
            else
            {
                cost.push_back(0);
            }
            ++num_vars;
        }
        // convert vector to array
        std::vector<size_t> coords(num_vars, 0);
        std::vector<size_t> shape(num_vars, (param_.max_number_objects + 1));
        //marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesP(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesA(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesD(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesT(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesDiv(shape.begin(), shape.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
            //std::cout << "...........................................StructuredLearningTrackingInferenceModel::add_detection_factors node state" << std::endl;
            if (param_.with_tracklets)
            {
                energy = 0;
                // add all detection factors of the internal nodes
                for (std::vector<Traxel>::const_iterator trax_it = tracklet_map_[n].begin();
                        trax_it != tracklet_map_[n].end(); ++trax_it)
                {
                    e = param_.detection(*trax_it, state);
                    energy += e;

                    energy += generateRandomOffset(Detection, e, *trax_it, 0, state);
                }
                // add all transition factors of the internal arcs
                Traxel tr_prev;
                bool first = true;
                for (std::vector<Traxel>::const_iterator trax_it = tracklet_map_[n].begin();
                        trax_it != tracklet_map_[n].end(); ++trax_it)
                {
                    LOG(logDEBUG4) << "internal arcs traxel " << *trax_it;
                    Traxel tr = *trax_it;
                    if (!first)
                    {
                        e = param_.transition( get_transition_probability(tr_prev, tr, state) );
                        energy += e;
                        energy += generateRandomOffset(Transition, e, tr_prev, tr);
                    }
                    else
                    {
                        first = false;
                    }
                    tr_prev = tr;
                }

            }
            else
            {
                e = param_.detection(traxel_map_[n], state);
                //std::cout << "detection=e   >" << e << std::endl;
                f = param_.detectionNoWeight(traxel_map_[n], state);
                //std::cout << "detection=e/w >" << f << std::endl;

                //w = param_.detection_weight;
                //std::cout << "detection=w   >" << w << std::endl;

                energy = f;

                energy += generateRandomOffset(Detection, e, traxel_map_[n], 0, state);
            }
            //std::cout << "===>" << energy << std::endl;
            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: detection[" << state
                           << "] = " << energy;
            for (size_t var_idx = 0; var_idx < num_vars; ++var_idx)
            {
                coords[var_idx] = state;
                // if only one of the variables is > 0, then it is an appearance in this time frame
                // or a disappearance in the next timeframe. Hence, add the cost of appearance/disappearance
                // to the detection cost


                //energies(coords.begin()) = conservationParam_.detection_weight * energy + state * cost[var_idx]; // state == m



                //std::cout << "___> " << conservationParam_.detection_weight << "   " << energy << "   " << state << "    " << cost[var_idx] << "   " << var_idx <<std::endl;
                energiesP(coords.begin()) = energy;




                if (var_idx == 0) // A
                    energiesA(coords.begin()) = state;// * cost[var_idx];
                else // D
                    energiesD(coords.begin()) = state;// * cost[var_idx];






                coords[var_idx] = 0;
                LOG(logDEBUG4) << "StructuredLearningTrackingInferenceModel::add_finite_factors: var_idx "
                               << var_idx << " = " << energy;
            }
            // also this energy if both variables have the same state
            if (num_vars == 2)
            {
                coords[0] = state;
                coords[1] = state;
                // only pay detection energy if both variables are on



                //energies(coords.begin()) = conservationParam_.detection_weight * energy;

                energiesP(coords.begin()) = energy;
                coords[0] = 0;
                coords[1] = 0;

                LOG(logDEBUG4) << "StructuredLearningTrackingInferenceModel::add_finite_factors: var_idxs 0 and var_idx 1 = "
                               << energy;
            }
        } // end for state

        //std::cout << "energies ==========================================================" << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energies(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << "energiesP ................." << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energiesP(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << "energiesA ................." << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energiesA(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << "energiesD ................." << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energiesD(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }

//        marray::Marray<double> errorMA(shape.begin(), shape.end(), param_.forbidden_cost);
//        errorMA = energies - conservationParam_.detection_weight * energiesP - cost[0]*energiesA - cost[1]*energiesD;



        // a decomposition of energies:
        // w_detection = conservationParam_.detection_weight
        // w_appearance = cost[0]
        // w_disappearance = cost[1]


//        std::cout << "errorMA ................." << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << errorMA(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }


//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                if (errorMA(coords.begin()) > 0.001 )
//                    std::cout << " We have a problem!!!";
//            }
//            //std::cout << std::endl;
//        }

        LOG(logDEBUG3) << "StructuredLearningTrackingInferenceModel::add_finite_factors: adding table to pgm";
        //functor add detection table
        //factorIndex = add_div_m_best_perturbation(energies, Detection, factorIndex);


        std::vector<size_t> weightIDs;
        weightIDs.push_back((size_t)0);
        weightIDs.push_back((size_t)3);
        weightIDs.push_back((size_t)4);
        std::vector<marray::Marray<double>> features;
        features.push_back(energiesP);
        features.push_back(energiesA);
        features.push_back(energiesD);

        std::vector<size_t> varShape;
        varShape.push_back((size_t)param_.max_number_objects+1);
        varShape.push_back((size_t)param_.max_number_objects+1);
        std::cout << varShape[0]*varShape[1] << " features[0].size()=" << features[0].size() << " numVar = " << num_vars << " maxnumobj+1= " << param_.max_number_objects + 1<< std::endl;
        opengm::functions::learnable::LSumOfExperts<double,size_t,size_t> funEnergies (varShape,weights_,weightIDs,features);

        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);

        sort(vi.begin(), vi.end());
        model_.addFactor(funcId, vi.begin(), vi.end());
        detection_f_node_map_[n] = model_.numberOfFactors() - 1;
    } // end for node n

    return factorIndex;
}

size_t StructuredLearningTrackingInferenceModel::add_transition_factors(const HypothesesGraph& g, size_t factorIndex)
{
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_transition_factors" << weights_.getWeight((size_t) 0) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_transition_factors" << weights_.getWeight((size_t) 1) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_transition_factors" << weights_.getWeight((size_t) 2) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_transition_factors" << weights_.getWeight((size_t) 3) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_transition_factors" << weights_.getWeight((size_t) 4) << std::endl;
    ////
    //// add transition factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());

    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_finite_factors: add transition factors";

    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        size_t vi[] = { arc_map_[a] };
        std::vector<size_t> coords(1, 0);

        std::vector<size_t> shape(1, (param_.max_number_objects + 1));
        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
            Traxel tr1, tr2;
            if (param_.with_tracklets)
            {
                tr1 = tracklet_map_[g.source(a)].back();
                tr2 = tracklet_map_[g.target(a)].front();
            }
            else
            {
                tr1 = traxel_map_[g.source(a)];
                tr2 = traxel_map_[g.target(a)];
            }

            double energy = param_.transition(get_transition_probability(tr1, tr2, state));
            energy += generateRandomOffset(Transition, energy, tr1, tr2);

            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: transition[" << state
                           << "] = " << energy;
            coords[0] = state;
            energies(coords.begin()) = energy;
            coords[0] = 0;
        }
        factorIndex = add_div_m_best_perturbation(energies, Transition, factorIndex);


        std::vector<size_t> weightIDs;
        weightIDs.push_back((size_t)2);
        std::vector<marray::Marray<double>> features;
        features.push_back(energies);

        std::vector<size_t> varShape;
        varShape.push_back((size_t)1+param_.max_number_objects);
        std::cout << varShape[0] << " features[0].size()=" << features[0].size() << " numVar = " << 1 << " maxnumobj+1= " << param_.max_number_objects + 1<< std::endl;
        opengm::functions::learnable::LSumOfExperts<double,size_t,size_t> funEnergies (varShape,weights_,weightIDs,features);

        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);
        model_.addFactor(funcId, vi, vi + 1);
    }

    return factorIndex;
}

size_t StructuredLearningTrackingInferenceModel::add_division_factors(const HypothesesGraph& g, size_t factorIndex)
{
    if(!param_.with_divisions)
    {
        return factorIndex;
    }

    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_division_factors" << weights_.getWeight((size_t) 0) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_division_factors" << weights_.getWeight((size_t) 1) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_division_factors" << weights_.getWeight((size_t) 2) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_division_factors" << weights_.getWeight((size_t) 3) << std::endl;
    std::cout << "...........................StructuredLearningTrackingInferenceModel::add_division_factors" << weights_.getWeight((size_t) 4) << std::endl;
    ////
    //// add division factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());

    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_finite_factors: add division factors";
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        if (div_node_map_.count(n) == 0)
        {
            continue;
        }
        size_t vi[] = { div_node_map_[n] };
        std::vector<size_t> coords(1, 0);
        std::vector<size_t> shape(1, 2);
        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= 1; ++state)
        {
            double energy;
            Traxel tr;
            if (param_.with_tracklets)
            {
                tr = tracklet_map_[n].back();
            }
            else
            {
                tr = traxel_map_[n];
            }
            energy = param_.division(tr, state);
            energy += generateRandomOffset(Division, energy,  tr);

            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: division[" << state
                           << "] = " << energy;
            coords[0] = state;
            energies(coords.begin()) = energy;
            coords[0] = 0;
        }
        factorIndex = add_div_m_best_perturbation(energies, Division, factorIndex);


        std::vector<size_t> weightIDs;
        weightIDs.push_back((size_t)4);
        std::vector<marray::Marray<double>> features;
        features.push_back(energies);

        std::vector<size_t> varShape;
        varShape.push_back((size_t)2);
        std::cout << varShape[0] << " features[0].size()=" << features[0].size() << " numVar = " << 1 << " maxnumobj+1= " << 1 + 1<< std::endl;
        opengm::functions::learnable::LSumOfExperts<double,size_t,size_t> funEnergies (varShape,weights_,weightIDs,features);

        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);
        model_.addFactor(funcId, vi, vi + 1);
    }

    return factorIndex;
}




} // namespace pgmlink

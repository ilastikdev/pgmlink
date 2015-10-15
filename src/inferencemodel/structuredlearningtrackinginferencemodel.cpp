#ifndef OPENGM_UNSIGNED_INTEGER_POW_HXX_
#define OPENGM_UNSIGNED_INTEGER_POW_HXX_
#endif

#include <opengm/functions/learnable/lweightedsum_of_functions.hxx>
#include "pgmlink/inferencemodel/structuredlearningtrackinginferencemodel.h"
#include <boost/python.hpp>
#include <limits>

namespace pgmlink
{

size_t StructuredLearningTrackingInferenceModel::add_detection_factors(const HypothesesGraph& g, size_t factorIndex)
{
    ////
    //// add detection factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());

    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_finite_factors: add detection factors";

    double minEnergyP = (double)std::numeric_limits<double>::infinity();
    double maxEnergyP = 0.0;
    double minEnergyA = (double)std::numeric_limits<double>::infinity();
    double maxEnergyA = 0.0;
    double minEnergyD = (double)std::numeric_limits<double>::infinity();
    double maxEnergyD = 0.0;

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        size_t num_vars = 0;
//        std::vector<size_t> vi;
//        std::vector<double> cost;

//        int node_begin_time = -1;
//        int node_end_time = -1;
//            node_begin_time = traxel_map_[n].Timestep;
//            node_end_time = traxel_map_[n].Timestep;

        double energy;//, e, f;//, w;
        if (app_node_map_.count(n) > 0)
        {
//            vi.push_back(app_node_map_[n]);
//            // "<" holds if there are only tracklets in the first frame
//            if (node_begin_time <= g.earliest_timestep())
//            {
//                // pay no appearance costs in the first timestep
//                cost.push_back(0.);
//            }
//            else
//            {
//                energy = param_.appearance_cost(traxel_map_[n]);
//                cost.push_back(energy);
//            }
            ++num_vars;
        }

        if (dis_node_map_.count(n) > 0)
        {
//            vi.push_back(dis_node_map_[n]);
//            // "<" holds if there are only tracklets in the last frame
//            if (node_end_time < g.latest_timestep())
//            {
//                energy = param_.disappearance_cost(traxel_map_[n]);

//                cost.push_back(energy);
//            }
//            else
//            {
//                cost.push_back(0);
//            }
            ++num_vars;
        }
        // convert vector to array
//        std::vector<size_t> coords(num_vars, 0);

//        std::vector<size_t> shape(num_vars, (param_.max_number_objects + 1));
//        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);
//        marray::Marray<double> energiesP(shape.begin(), shape.end(), param_.forbidden_cost);
//        marray::Marray<double> energiesA(shape.begin(), shape.end(), param_.forbidden_cost);
//        marray::Marray<double> energiesD(shape.begin(), shape.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
                //e = param_.detection(traxel_map_[n], state);
                //f = param_.detectionNoWeight(traxel_map_[n], state);
                energy = param_.detectionNoWeight(traxel_map_[n], state);

            for (size_t var_idx = 0; var_idx < num_vars; ++var_idx)
            {
//                coords[var_idx] = state;

//                energies(coords.begin()) = inferenceWeights_.getWeight((size_t)0) * energy + state * cost[var_idx]; // state == m
//                energiesP(coords.begin()) = energy;

                if (energy < minEnergyP)
                    minEnergyP = energy;
                if (energy > maxEnergyP)
                    maxEnergyP = energy;

                if (var_idx == 0){ // A
//                    energiesA(coords.begin()) = state;// * cost[var_idx];
//                    energiesAVec(coordsVec.begin()) = state;// * cost[var_idx];
                    if (state < minEnergyA)
                        minEnergyA = state;
                    if (state > maxEnergyA)
                        maxEnergyA = state;
                }
                else{ // D
//                    energiesD(coords.begin()) = state;// * cost[var_idx];
//                    energiesDVec(coordsVec.begin()) = state;// * cost[var_idx];
                    if (state < minEnergyD)
                        minEnergyD = state;
                    if (state > maxEnergyD)
                        maxEnergyD = state;
                }

//                coords[var_idx] = 0;
            }
            // also this energy if both variables have the same state
            if (num_vars == 2)
            {
//                coords[0] = state;
//                coords[1] = state;

//                energiesP(coords.begin()) = energy;

                if (energy < minEnergyP)
                    minEnergyP = energy;
                if (energy > maxEnergyP)
                    maxEnergyP = energy;

//                coords[0] = 0;
//                coords[1] = 0;

            }
        } // end for state
    } // end for node n

    std::cout << "minEnergyDet " << minEnergyP << std::endl;
    std::cout << "maxEnergyDet " << maxEnergyP << std::endl;
    std::cout << "minEnergyApp " << minEnergyA << std::endl;
    std::cout << "maxEnergyApp " << maxEnergyA << std::endl;
    std::cout << "minEnergyDis " << minEnergyD << std::endl;
    std::cout << "maxEnergyDis " << maxEnergyD << std::endl;

    double minEnergyPnew = (double)std::numeric_limits<double>::infinity();
    double maxEnergyPnew = 0.0;
    double minEnergyAnew = (double)std::numeric_limits<double>::infinity();
    double maxEnergyAnew = 0.0;
    double minEnergyDnew = (double)std::numeric_limits<double>::infinity();
    double maxEnergyDnew = 0.0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        size_t num_vars = 0;
        std::vector<size_t> vi;
        std::vector<double> cost;

        int node_begin_time = -1;
        int node_end_time = -1;
            node_begin_time = traxel_map_[n].Timestep;
            node_end_time = traxel_map_[n].Timestep;

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
                energy = param_.appearance_cost(traxel_map_[n]);
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
                energy = param_.disappearance_cost(traxel_map_[n]);

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
//        std::vector<size_t> coordsSym(num_vars, 0);

        std::vector<size_t> shape(num_vars, (param_.max_number_objects + 1));
        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesP(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesA(shape.begin(), shape.end(), param_.forbidden_cost);
        marray::Marray<double> energiesD(shape.begin(), shape.end(), param_.forbidden_cost);

//        std::vector<size_t> coordsVec(1, 0);
//        std::vector<size_t> shapeVec(1,(param_.max_number_objects + 1)*(param_.max_number_objects + 1));
//        marray::Marray<double> energiesVec(shapeVec.begin(), shapeVec.end(), param_.forbidden_cost);
//        marray::Marray<double> energiesPVec(shapeVec.begin(), shapeVec.end(), param_.forbidden_cost);
//        marray::Marray<double> energiesAVec(shapeVec.begin(), shapeVec.end(), param_.forbidden_cost);
//        marray::Marray<double> energiesDVec(shapeVec.begin(), shapeVec.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
                e = param_.detection(traxel_map_[n], state);
                f = param_.detectionNoWeight(traxel_map_[n], state);
                energy = f;

            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: detection[" << state
                           << "] = " << energy;
            for (size_t var_idx = 0; var_idx < num_vars; ++var_idx)
            {
                coords[var_idx] = state;
//                coordsVec[0] = coords[0]*(param_.max_number_objects + 1)+coords[1];


                // if only one of the variables is > 0, then it is an appearance in this time frame
                // or a disappearance in the next timeframe. Hence, add the cost of appearance/disappearance
                // to the detection cost


                energies(coords.begin()) = inferenceWeights_.getWeight((size_t)0) * energy + state * cost[var_idx]; // state == m
                energiesP(coords.begin()) = (energy-minEnergyP)/(maxEnergyP-minEnergyP);

//                energiesVec(coordsVec.begin()) = inferenceWeights_.getWeight((size_t)0) * energy + state * cost[var_idx]; // state == m
//                energiesPVec(coordsVec.begin()) = energy;

                if (energiesP(coords.begin()) < minEnergyPnew)
                    minEnergyPnew = energiesP(coords.begin());
                if (energiesP(coords.begin()) > maxEnergyPnew)
                    maxEnergyPnew = energiesP(coords.begin());

                if (var_idx == 0){ // A
                    energiesA(coords.begin()) = (state-minEnergyA)/(maxEnergyA-minEnergyA);// * cost[var_idx];
//                    energiesAVec(coordsVec.begin()) = state;// * cost[var_idx];

                    if (energiesA(coords.begin()) < minEnergyAnew)
                        minEnergyAnew = energiesA(coords.begin());
                    if (energiesA(coords.begin()) > maxEnergyAnew)
                        maxEnergyAnew = energiesA(coords.begin());
                }
                else{ // D
                    energiesD(coords.begin()) = (state-minEnergyD)/(maxEnergyD-minEnergyD);// * cost[var_idx];
//                    energiesDVec(coordsVec.begin()) = state;// * cost[var_idx];

                    if (energiesD(coords.begin()) < minEnergyDnew)
                        minEnergyDnew = energiesD(coords.begin());
                    if (energiesD(coords.begin()) > maxEnergyDnew)
                        maxEnergyDnew = energiesD(coords.begin());
                }

                coords[var_idx] = 0;
                LOG(logDEBUG4) << "StructuredLearningTrackingInferenceModel::add_finite_factors: var_idx "
                               << var_idx << " = " << energy;
            }
            // also this energy if both variables have the same state
            if (num_vars == 2)
            {
                coords[0] = state;
                coords[1] = state;
//                coordsVec[0] = coords[0]*(param_.max_number_objects + 1)+coords[1];
                // only pay detection energy if both variables are on

                energies(coords.begin()) = inferenceWeights_.getWeight((size_t)0) * energy;
                energiesP(coords.begin()) = (energy-minEnergyP)/(maxEnergyP-minEnergyP);

//                energiesVec(coordsVec.begin()) = inferenceWeights_.getWeight((size_t)0) * energy;
//                energiesPVec(coordsVec.begin()) = energy;

                if (energiesP(coords.begin()) < minEnergyPnew)
                    minEnergyPnew = energiesP(coords.begin());
                if (energiesP(coords.begin()) > maxEnergyPnew)
                    maxEnergyPnew = energiesP(coords.begin());

                coords[0] = 0;
                coords[1] = 0;

                LOG(logDEBUG4) << "StructuredLearningTrackingInferenceModel::add_finite_factors: var_idxs 0 and var_idx 1 = "
                               << energy;
            }
        } // end for state

//        marray::Marray<double> errorMA(shape.begin(), shape.end(), param_.forbidden_cost);
//        errorMA = energies - inferenceWeights_.getWeight((size_t)0) * energiesP - cost[0]*energiesA - cost[1]*energiesD;

//        marray::Marray<double> errorMAVec(shapeVec.begin(), shapeVec.end(), param_.forbidden_cost);
//        errorMAVec = energiesVec - inferenceWeights_.getWeight((size_t)0) * energiesPVec - cost[0]*energiesAVec - cost[1]*energiesDVec;

//        std::cout << "energies ==========================================================" << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energies(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }
//        for (size_t i = 0; i < (param_.max_number_objects+1)*(param_.max_number_objects+1); ++i){
//            coordsVec[0] = i;
//            std::cout << energiesVec(coordsVec.begin()) << " ";
//        }
//        std::cout << std::endl;
//        std::cout << "energiesDet .........DETECTION........" << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                //if(coords[1]>0) energiesP(coords.begin()) = 0;
//                //if(coords[0]==0) energiesP(coords.begin()) = 0;
//                std::cout <<  energiesP(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }
//        for (size_t i = 0; i < (param_.max_number_objects+1)*(param_.max_number_objects+1); ++i){
//            coordsVec[0] = i;
//            //if(coords[1]>3) energiesPVec(coordsVec.begin()) = 0;
//            std::cout << energiesPVec(coordsVec.begin()) << " ";
//        }
//        std::cout << std::endl;
//        std::cout << "energiesApp ..........APPEARANCE......." << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energiesA(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }
//        for (size_t i = 0; i < (param_.max_number_objects+1)*(param_.max_number_objects+1); ++i){
//            coordsVec[0] = i;
//            std::cout << energiesAVec(coordsVec.begin()) << " ";
//        }
//        std::cout << std::endl;
//        std::cout << "energiesDis .........DISAPPEARANCE........" << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                std::cout << energiesD(coords.begin()) << " ";
//            }
//            std::cout << std::endl;
//        }
//        for (size_t i = 0; i < (param_.max_number_objects+1)*(param_.max_number_objects+1); ++i){
//            coordsVec[0] = i;
//            std::cout << energiesDVec(coordsVec.begin()) << " ";
//        }
//        std::cout << "END print" << std::endl;
//        std::cout << "errorMA ................." << std::endl;
//        bool noErrorFlag = true;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            for (size_t state2 = 0; state2 <= param_.max_number_objects; ++state2){
//                coords[0] = state;
//                coords[1] = state2;
//                //std::cout << errorMA(coords.begin()) << " ";
//                if(errorMA(coords.begin()) >  0.00000000000001)
//                    noErrorFlag &= false;
//            }
//            //std::cout << std::endl;
//        }
//        if(!noErrorFlag)
//            std::cout << "==================================================================>ERROR in MA" << std::endl;

//        std::cout << "param_.detection_weight = " << param_.detection_weight << std::endl;
//        for (int i = 0; i < 5; ++i)
//            std::cout << "===== inf weights ===> " << inferenceWeights_.getWeight((size_t)i) << std::endl;

        // a decomposition of energies:
        // w_detection = conservationParam_.detection_weight
        // w_appearance = cost[0]
        // w_disappearance = cost[1]

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
//        std::cout << "funEnergies" << std::endl;
        opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t> funEnergies (varShape,inferenceWeights_,weightIDs,features);

//        std::vector<marray::Marray<double>> featuresVec;
//        featuresVec.push_back(energiesPVec);
//        featuresVec.push_back(energiesAVec);
//        featuresVec.push_back(energiesDVec);
//        std::vector<size_t> varShapeVec;
//        varShapeVec.push_back((size_t)(param_.max_number_objects+1)*(param_.max_number_objects+1));
//        std::cout << "--->funEnergies" << std::endl;
//        opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t> funEnergies (varShapeVec,inferenceWeights_,weightIDs,featuresVec);

//        std::cout << "addFunction" << std::endl;
        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);

//        for (auto i = vi.begin(); i !=vi.end(); ++i){
//            std::cout << "i=" << *i << std::endl;
//        }

//        std::cout << "sort" << std::endl;
        sort(vi.begin(), vi.end());
//        std::cout << "addFactor" << std::endl;
        model_.addFactor(funcId, vi.begin(), vi.end());
//        std::cout << "numberOfFactors" << std::endl;
        detection_f_node_map_[n] = model_.numberOfFactors() - 1;
//        std::cout << "end node" << std::endl;

    } // end for node n

    std::cout << "minEnergyDet new " << minEnergyPnew << std::endl;
    std::cout << "maxEnergyDet new " << maxEnergyPnew << std::endl;
    std::cout << "minEnergyApp new " << minEnergyAnew << std::endl;
    std::cout << "maxEnergyApp new " << maxEnergyAnew << std::endl;
    std::cout << "minEnergyDis new " << minEnergyDnew << std::endl;
    std::cout << "maxEnergyDis new " << maxEnergyDnew << std::endl;
    return factorIndex;
}

size_t StructuredLearningTrackingInferenceModel::add_transition_factors(const HypothesesGraph& g, size_t factorIndex)
{
    ////
    //// add transition factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());

    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_finite_factors: add transition factors";

    double minEnergy = (double)std::numeric_limits<double>::infinity();
    double maxEnergy = 0.0;

    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
//        size_t vi[] = { arc_map_[a] };
//        std::vector<size_t> coords(1, 0);

//        std::vector<size_t> shape(1, (param_.max_number_objects + 1));
//        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
            Traxel tr1, tr2;
            tr1 = traxel_map_[g.source(a)];
            tr2 = traxel_map_[g.target(a)];

            //double energy = param_.transition(get_transition_probability(tr1, tr2, state));
            double energy = param_.transition(tr1, tr2, state);
            //energy += generateRandomOffset(Transition, energy, tr1, tr2);

            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: transition[" << state
                           << "] = " << energy;
//            coords[0] = state;
//            energies(coords.begin()) = energy;
            if (energy < minEnergy)
                minEnergy = energy;
            if (energy > maxEnergy)
                maxEnergy = energy;
//            coords[0] = 0;
        }
        //factorIndex = add_div_m_best_perturbation(energies, Transition, factorIndex);


//        std::vector<size_t> weightIDs;
//        weightIDs.push_back((size_t)2);
//        std::vector<marray::Marray<double>> features;
//        features.push_back(energies);

//        std::vector<size_t> varShape;
//        varShape.push_back((size_t)1+param_.max_number_objects);

//        opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t> funEnergies (varShape,inferenceWeights_,weightIDs,features);

//        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);
//        model_.addFactor(funcId, vi, vi + 1);

//        std::cout << "energies ==========================TRANSITION================================" << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            coords[0] = state;
//            std::cout << energies(coords.begin()) << " ";
//        }
//        std::cout << std::endl;
    }

    std::cout << "minEnergyTran  " << minEnergy << std::endl;
    std::cout << "maxEnergyTran " << maxEnergy << std::endl;

    double minEnergyNew = (double)std::numeric_limits<double>::infinity();
    double maxEnergyNew = 0.0;

    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        size_t vi[] = { arc_map_[a] };
        std::vector<size_t> coords(1, 0);

        std::vector<size_t> shape(1, (param_.max_number_objects + 1));
        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= param_.max_number_objects; ++state)
        {
            Traxel tr1, tr2;
                tr1 = traxel_map_[g.source(a)];
                tr2 = traxel_map_[g.target(a)];

            //double energy = param_.transition(get_transition_probability(tr1, tr2, state));
            double energy = param_.transition(tr1, tr2, state);
            //energy += generateRandomOffset(Transition, energy, tr1, tr2);

            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: transition[" << state
                           << "] = " << energy;
            coords[0] = state;
            energies(coords.begin()) = (energy-minEnergy)/(maxEnergy-minEnergy);

            if (energies(coords.begin()) < minEnergyNew)
                minEnergyNew = energies(coords.begin());
            if (energies(coords.begin()) > maxEnergyNew)
                maxEnergyNew = energies(coords.begin());
            coords[0] = 0;
        }
        //factorIndex = add_div_m_best_perturbation(energies, Transition, factorIndex);


        std::vector<size_t> weightIDs;
        weightIDs.push_back((size_t)2);
        std::vector<marray::Marray<double>> features;
        features.push_back(energies);

        std::vector<size_t> varShape;
        varShape.push_back((size_t)1+param_.max_number_objects);

        opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t> funEnergies (varShape,inferenceWeights_,weightIDs,features);

        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);
        model_.addFactor(funcId, vi, vi + 1);

//        std::cout << "energies ==========================TRANSITION================================" << std::endl;
//        for (size_t state = 0; state <= param_.max_number_objects; ++state){
//            coords[0] = state;
//            std::cout << energies(coords.begin()) << " ";
//        }
//        std::cout << std::endl;
    }

    std::cout << "minEnergyTran new " << minEnergyNew << std::endl;
    std::cout << "maxEnergyTran new " << maxEnergyNew << std::endl;
    return factorIndex;
}

size_t StructuredLearningTrackingInferenceModel::add_division_factors(const HypothesesGraph& g, size_t factorIndex)
{
    if(!param_.with_divisions)
    {
        return factorIndex;
    }
    ////
    //// add division factors
    ////
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map_ = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map_ =
        g.get(node_tracklet());

    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_finite_factors: add division factors";

    double minEnergy = (double)std::numeric_limits<double>::infinity();
    double maxEnergy = 0.0;

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        if (div_node_map_.count(n) == 0)
        {
            continue;
        }
//        size_t vi[] = { div_node_map_[n] };
//        std::vector<size_t> coords(1, 0);
//        std::vector<size_t> shape(1, 2);
//        marray::Marray<double> energies(shape.begin(), shape.end(), param_.forbidden_cost);

        for (size_t state = 0; state <= 1; ++state)
        {
            double energy;
            Traxel tr;
            tr = traxel_map_[n];
            //energy = param_.division(tr, state);
            energy = param_.divisionNoWeight(tr, state);


            //energy += generateRandomOffset(Division, energy,  tr);

//            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: division[" << state
//                           << "] = " << energy;
//            std::cout << "StructuredLearningTrackingInferenceModel::add_finite_factors: tr=" << tr << " state=" << state << " division[" << state
//                           << "] = " << energy << std::endl;
//            coords[0] = state;
//            energies(coords.begin()) = energy;

            if (energy < minEnergy)
                minEnergy = energy;
            if (energy > maxEnergy)
                maxEnergy = energy;
//            coords[0] = 0;
        }
        //factorIndex = add_div_m_best_perturbation(energies, Division, factorIndex);

//        std::vector<size_t> weightIDs;
//        weightIDs.push_back((size_t)1);
//        std::vector<marray::Marray<double>> features;
//        features.push_back(energies);
//        std::vector<size_t> varShape;
//        varShape.push_back((size_t)2);

//        opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t> funEnergies (varShape,inferenceWeights_,weightIDs,features);

//        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);
//        model_.addFactor(funcId, vi, vi + 1);

//        std::cout << "=========================DIVISION=================================" << std::endl;
//        for (size_t state = 0; state < 2; ++state){
//            coords[0] = state;
//            std::cout << energies(coords.begin()) << " ";
//        }
//        std::cout << std::endl;
    }

    std::cout << "minEnergyDiv  " << minEnergy << std::endl;
    std::cout << "maxEnergyDiv " << maxEnergy << std::endl;

    double minEnergyNew = (double)std::numeric_limits<double>::infinity();
    double maxEnergyNew = 0.0;

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
            tr = traxel_map_[n];
            //energy = param_.division(tr, state);
            energy = param_.divisionNoWeight(tr, state);


            //energy += generateRandomOffset(Division, energy,  tr);

            LOG(logDEBUG2) << "StructuredLearningTrackingInferenceModel::add_finite_factors: division[" << state
                           << "] = " << energy;
//            std::cout << "StructuredLearningTrackingInferenceModel::add_finite_factors: tr=" << tr << " state=" << state << " division[" << state
//                           << "] = " << energy << std::endl;
            coords[0] = state;
            energies(coords.begin()) = (energy-minEnergy)/(maxEnergy-minEnergy);

            if (energies(coords.begin()) < minEnergyNew)
                minEnergyNew = energies(coords.begin());
            if (energies(coords.begin()) > maxEnergyNew)
                maxEnergyNew = energies(coords.begin());
            coords[0] = 0;
        }
        //factorIndex = add_div_m_best_perturbation(energies, Division, factorIndex);

        std::vector<size_t> weightIDs;
        weightIDs.push_back((size_t)1);
        std::vector<marray::Marray<double>> features;
        features.push_back(energies);
        std::vector<size_t> varShape;
        varShape.push_back((size_t)2);

        opengm::functions::learnable::LWeightedSumOfFunctions<double,size_t,size_t> funEnergies (varShape,inferenceWeights_,weightIDs,features);

        typename GraphicalModelType::FunctionIdentifier funcId = model_.addFunction(funEnergies);
        model_.addFactor(funcId, vi, vi + 1);

//        std::cout << "=========================DIVISION=================================" << std::endl;
//        for (size_t state = 0; state < 2; ++state){
//            coords[0] = state;
//            std::cout << energies(coords.begin()) << " ";
//        }
//        std::cout << std::endl;
    }

    std::cout << "minEnergyDiv New " << minEnergyNew << std::endl;
    std::cout << "maxEnergyDiv New " << maxEnergyNew << std::endl;

    return factorIndex;
}
} // namespace pgmlink

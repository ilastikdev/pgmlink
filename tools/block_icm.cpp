
#include <iostream>
#include <algorithm>

#include "pgmlink/pgm.h"
#include "pgmlink/constraint_pool.hxx"
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>
#include <opengm/inference/icm.hxx>

std::vector<size_t> map_factor_indices(std::map<size_t, size_t>& index_mapping, pgmlink::pgm::OpengmModelDeprecated::ogmGraphicalModel& model, size_t factor_id)
{
    std::vector<size_t> new_indices;

    for(size_t i = 0; i < model[factor_id].numberOfVariables(); i++)
    {
        if(index_mapping.find(model[factor_id].variableIndex(i)) == index_mapping.end())
            return std::vector<size_t>();
        else
            new_indices.push_back(index_mapping[model[factor_id].variableIndex(i)]);
    }

    return new_indices;
}

int main(int argc, char** argv)
{
    if(argc != 3)
    {
        std::cout << "Perform something similar to block-ICM on stored conservation tracking models with constraints. 2014 (c) Carsten Haubold" << std::endl;
        std::cout << "\nUSAGE: " << argv[0] << " model.h5 constraints.cp" << std::endl;
        return 0;
    }

    std::string filename_model(argv[1]);
    std::string filename_constraints(argv[2]);

    // load model and constraints from disk
    pgmlink::pgm::OpengmModelDeprecated::ogmGraphicalModel model;
    opengm::hdf5::load(model, filename_model, "model");

    std::ifstream constraint_pool_input(filename_constraints);
    std::map< size_t, std::vector<size_t> > nodes_per_timestep;
    pgmlink::pgm::ConstraintPool cp;

    {
        boost::archive::text_iarchive ia(constraint_pool_input);
        ia & cp;
        ia & nodes_per_timestep;
    }
    constraint_pool_input.close();

    // dump statistics
    std::cout << "Loaded Model from " << filename_model << std::endl;
    std::cout << "\tVariables: " << model.numberOfVariables() << std::endl;
    std::cout << "\tFactors: " << model.numberOfFactors() << std::endl;

    std::cout << "\nLoaded ConstraintPool from " << filename_constraints << std::endl;
    std::cout << "\tNum Constraints: " << cp.get_num_constraints() << std::endl;

    std::vector<size_t> solution(model.numberOfVariables(), 0);

    // iterate over time steps, create submodel (view?) just for pairwise time frames, CPLEX inference on submodels only
    for(auto timestep_it = nodes_per_timestep.begin(); timestep_it != nodes_per_timestep.end(); ++timestep_it)
    {
        std::vector<size_t>& nodes = timestep_it->second;
        // make sure all index relations stay
        std::sort(nodes.begin(), nodes.end(), std::greater<size_t>());

        std::map<size_t, size_t> index_mapping;
        pgmlink::pgm::OpengmModelDeprecated::ogmGraphicalModel submodel;

        // copy all nodes and store a mapping
        for(auto node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
        {
            index_mapping[*node_it] = submodel.addVariable(model.numberOfLabels(*node_it));
        }

        // copy all contained factors
        for(size_t factor_id = 0; factor_id < model.numberOfFactors(); ++factor_id)
        {
            std::vector<size_t> remapped_factor_indices = map_factor_indices(index_mapping, model, factor_id);
            if(model[factor_id].numberOfVariables() == 0 || remapped_factor_indices.size() > 0)
            {
                pgmlink::pgm::OpengmModelDeprecated::ogmGraphicalModel::IndependentFactorType independent_factor(model[factor_id]);
                pgmlink::pgm::OpengmModelDeprecated::FunctionIdentifier sub_id = submodel.addFunction(independent_factor.function());
                submodel.addFactor(sub_id, remapped_factor_indices.begin(), remapped_factor_indices.end());
            }
        }

        // create inference type and add constraints
        opengm::LPCplex<pgmlink::pgm::OpengmModelDeprecated::ogmGraphicalModel, pgmlink::pgm::OpengmModelDeprecated::ogmAccumulator>::Parameter param;
        param.verbose_ = true;
        param.integerConstraint_ = true;
        param.epGap_ = 0.0;
        opengm::LPCplex<pgmlink::pgm::OpengmModelDeprecated::ogmGraphicalModel, pgmlink::pgm::OpengmModelDeprecated::ogmAccumulator> inf(submodel, param);
        cp.add_constraints_to_problem(submodel, inf, index_mapping);

        // infer
        opengm::InferenceTermination status = inf.infer();

        if (status != opengm::NORMAL) {
            throw std::runtime_error("CPLEX optimizer terminated abnormally");
        }

        // extract and print solution
        std::vector<size_t> subsolution;
        status = inf.arg(solution);

        if (status != opengm::NORMAL) {
            throw std::runtime_error("Could not extract solution from CPLEX");
        }

        // map solution back
        for(auto node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
        {
            solution[*node_it] = subsolution[index_mapping[*node_it]];
        }

        // evaluate model with current solution
        std::cout << "Solution after optimizing timestep " << timestep_it->first << " has energy: " << model.evaluate(solution);
    }

    std::cout << "Done";
}

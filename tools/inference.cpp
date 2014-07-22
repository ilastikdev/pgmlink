#include <iostream>
#include "pgmlink/pgm.h"
#include "pgmlink/constraint_pool.hxx"
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>
#include <opengm/inference/icm.hxx>

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		std::cout << "Inference runner on stored conservation tracking models with constraints. 2014 (c) Carsten Haubold" << std::endl;
		std::cout << "\nUSAGE: " << argv[0] << " model.h5 constraints.cp CPLEX|ICM" << std::endl;
		return 0;
	}

	std::string filename_model(argv[1]);
	std::string filename_constraints(argv[2]);
	bool use_cplex = std::string(argv[3]) == "CPLEX";

	// load model and constraints from disk
	pgmlink::pgm::OpengmModelDeprecated::ogmGraphicalModel model;
    opengm::hdf5::load(model, filename_model, "model");
    
    std::ifstream constraint_pool_input(filename_constraints);
    pgmlink::pgm::ConstraintPool cp;
    
    {
        boost::archive::text_iarchive ia(constraint_pool_input);
        ia & cp;
    }
    constraint_pool_input.close();

    // dump statistics
    std::cout << "Loaded Model from " << filename_model << std::endl;
    std::cout << "\tVariables: " << model.numberOfVariables() << std::endl;
    std::cout << "\tFactors: " << model.numberOfFactors() << std::endl;
    
    std::cout << "\nLoaded ConstraintPool from " << filename_constraints << std::endl;
    std::cout << "\tNum Constraints: " << cp.get_num_constraints() << std::endl;

    std::vector<pgmlink::pgm::OpengmModelDeprecated::ogmInference::LabelType> solution;

    if(use_cplex)
    {
    	opengm::LPCplex<pgmlink::pgm::OpengmModelDeprecated::ogmGraphicalModel, pgmlink::pgm::OpengmModelDeprecated::ogmAccumulator>::Parameter param;
	    param.verbose_ = true;
	    param.integerConstraint_ = true;
	    param.epGap_ = 0.0;
    	opengm::LPCplex<pgmlink::pgm::OpengmModelDeprecated::ogmGraphicalModel, pgmlink::pgm::OpengmModelDeprecated::ogmAccumulator> inf(model, param);
    	cp.add_constraints_to_problem(model, inf);

    	opengm::InferenceTermination status = inf.infer();

    	if (status != opengm::NORMAL) {
	        throw std::runtime_error("CPLEX optimizer terminated abnormally");
	    }

	    // extract and print solution
    	status = inf.arg(solution);

    	if (status != opengm::NORMAL) {
	        throw std::runtime_error("Could not extract solution from CPLEX");
	    }
    }
    else
    {
    	// configure ICM
	    opengm::ICM<pgmlink::pgm::OpengmModelDeprecated::ogmGraphicalModel, pgmlink::pgm::OpengmModelDeprecated::ogmAccumulator> inf(model);
	    cp.add_constraints_to_problem(model, inf);
	    cp.set_big_m(10000000.0);

	    // run inference
	    opengm::InferenceTermination status = inf.infer();

	    if (status != opengm::NORMAL) {
	        throw std::runtime_error("ICM optimizer terminated abnormally");
	    }

	    // extract and print solution
    	status = inf.arg(solution);

    	if (status != opengm::NORMAL) {
	        throw std::runtime_error("Could not extract solution from ICM");
	    }

    }
    std::cout << "\n===============================\nFound Solution:" << std::endl;

    for(auto it = solution.begin(); it != solution.end(); ++it)
    {
    	std::cout << *it << " ";
    }
    std::cout << std::endl;
}
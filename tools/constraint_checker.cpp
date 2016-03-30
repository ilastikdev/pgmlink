
#include <iostream>
#include <algorithm>
#include <string>

#include <boost/algorithm/string.hpp>
#include "pgmlink/pgm.h"
#include "pgmlink/inferencemodel/constraint_pool.hxx"
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>
#include <opengm/inference/icm.hxx>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

// #include <vigra/timing.hxx>

typedef pgmlink::PertGmType GraphicalModel;
typedef opengm::LPCplex<GraphicalModel, pgmlink::pgm::OpengmModelDeprecated::ogmAccumulator> CPLEX;
typedef opengm::LPCplex2<GraphicalModel, pgmlink::pgm::OpengmModelDeprecated::ogmAccumulator> CPLEX2;

typedef std::map<size_t, size_t> IndexMapping;
typedef std::vector<size_t> Solution;

int main(int argc, char** argv)
{
    // USETICTOC
    
    if(argc != 4)
    {
        std::cout << "Load a model, its constraints and a solution and check which constraints are violated. 2015 (c) Carsten Haubold" << std::endl;
        std::cout << "\nUSAGE: " << argv[0] << " model.h5 constraints.cp solution.txt(or solution as string)" << std::endl;
        return 0;
    }

    // TIC

    std::string filename_model(argv[1]);
    std::string filename_constraints(argv[2]);
    std::string filename_solution(argv[3]);

    double big_m = 10000000.0;

    // load model and constraints from disk
    GraphicalModel model;
    opengm::hdf5::load(model, filename_model, "conservationTracking");

    std::ifstream constraint_pool_input(filename_constraints);
    if(!constraint_pool_input.good())
        throw std::runtime_error("Couldn't open constraint file");

    std::map< size_t, std::vector<size_t> > nodes_per_timestep;
    pgmlink::pgm::ConstraintPool cp;
    {
        boost::archive::text_iarchive ia(constraint_pool_input);
        ia & cp;
        ia & nodes_per_timestep;
    }
    constraint_pool_input.close();
    cp.set_big_m(big_m);

    // dump statistics
    std::cout << "Loaded Model from " << filename_model << std::endl;
    std::cout << "\tVariables: " << model.numberOfVariables() << std::endl;
    std::cout << "\tFactors: " << model.numberOfFactors() << std::endl;

    std::cout << "Loaded ConstraintPool from " << filename_constraints << std::endl;
    std::cout << "\tNum Constraints: " << cp.get_num_constraints() << std::endl;

    Solution solution;
    std::ifstream solution_input(filename_solution);
    if(solution_input.good())
    {
        {
            boost::archive::text_iarchive ia(solution_input);
            ia & solution;
        }
        solution_input.close();
        std::cout << "Loaded solution from " << filename_solution << std::endl;
    }
    else
    {
        // assume input is a string:
        std::vector<std::string> labels;
        boost::split(labels, filename_solution, boost::is_any_of(", "), boost::token_compress_on);
        solution.resize(labels.size());
        std::transform(labels.begin(), labels.end(), solution.begin(), [](const std::string& s){
            return std::stol(s);
        });
        std::cout << "Loaded solution from input string" << std::endl;
    }
    std::cout << "\tNum variable labels: " << solution.size() << std::endl;

    assert(solution.size() == model.numberOfVariables());

    // evaluate
    double pure_value = model.evaluate(solution);
    

    opengm::ICM<GraphicalModel, pgmlink::pgm::OpengmModelDeprecated::ogmAccumulator> inf(model);
    cp.add_constraints_to_problem(model, inf);

    double value_with_big_m = model.evaluate(solution);

    std::cout << "\n--------------------------------\n" << std::endl;
    std::cout << "Model did violate: " << (value_with_big_m - pure_value) / big_m << " constraints" << std::endl;

}

/**
   @file
   @ingroup tracking
   @brief field-of-view filtering of a traxelstore
*/

#ifndef UNCERTAINTYPARAMETER_H
#define UNCERTAINTYPARAMETER_H
#include <iostream>
#include "pgmlink/pgmlink_export.h"
using namespace std;

namespace pgmlink {
	enum PGMLINK_EXPORT DistrId {
		GaussianPertubation,PerturbAndMAP,DiverseMbest,MbestCPLEX,ClassifierUncertainty
	};

	class PGMLINK_EXPORT UncertaintyParameter{
	public:
		std::size_t numberOfIterations;
		DistrId distributionId;
		std::vector<double> distributionParam;

		UncertaintyParameter(){
			numberOfIterations=1;

			distributionId=GaussianPertubation;
			//distributionId table:
			//0: Gauss normal
			//1: Gumbel Perturb&MAP
			//2: diverse-m-best deterministic
			//3: diverse-m-best cplex

			distributionParam = std::vector<double>(1,0);
			//various distribution parameters:
			//for distributions depending on a parameter
			//distributionParam holds a parameter value for
			//perturbing each of the following independently
			// - appearance
			// - disappearance
			// - detection
			// - transition
			// - division
		}
		UncertaintyParameter(std::size_t nOI,DistrId dI,std::vector<double> dP){
			numberOfIterations=nOI;
			distributionId = dI;
			distributionParam = dP;
		}
		//constructor for single-parameter distributions
		UncertaintyParameter(std::size_t nOI,DistrId dI,double dP){
			numberOfIterations=nOI;
			distributionId = dI;
			distributionParam = std::vector<double>(1,dP);

		}
		void print(){
			cout<<"number of iterations "<<numberOfIterations<<endl;
			cout<<"distribution Id "<<distributionId<<endl;
			cout<<"distribution Parameters "<<endl;
			for(std::vector<double>::iterator it = distributionParam.begin(); it != distributionParam.end(); ++it) {
				cout << *it << ", ";
			}
			cout<<endl;
		}
	};

} /* namespace pgmlink */

#endif /* UNCERTAINTYPARAMETER_H */

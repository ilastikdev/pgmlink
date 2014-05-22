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
		GaussianPertubation,PerturbAndMAP,DiverseMbest,MbestCPLEX
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
			cout<<"number of iterations "<<numberOfIterations;
			cout<<"distribution Id "<<distributionId;
			cout<<"distribution Parameters ";
			for(std::vector<double>::iterator it = distributionParam.begin(); it != distributionParam.end(); ++it) {
				cout << *it << ", ";
			}
		}
	};

} /* namespace pgmlink */

#endif /* UNCERTAINTYPARAMETER_H */

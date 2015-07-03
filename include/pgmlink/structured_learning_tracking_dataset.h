#ifndef OPENGM_DATASET_H
#define OPENGM_DATASET_H

#include <opengm/learning/dataset/dataset.hxx>
#include "pgmlink/inferencemodel/constrackinginferencemodel.h"
#include "pgmlink/tracking.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/graph.h"
#include "pgmlink/field_of_view.h"

namespace opengm {
namespace datasets{

template<class GM, class LOSS>
class StructuredLearningTrackingDataset : public Dataset<GM,LOSS,GM>{
public:
   typedef GM                     GMType;
   typedef GM                     GMWITHLOSS;
   typedef LOSS                   LossType;
   typedef typename GM::ValueType ValueType;
   typedef typename GM::IndexType IndexType;
   typedef typename GM::LabelType LabelType;
   typedef opengm::learning::Weights<ValueType> Weights;

   StructuredLearningTrackingDataset()
   {}

   StructuredLearningTrackingDataset(
       pgmlink::ConsTrackingInferenceModel::IndexType numModels,
       std::vector<pgmlink::FieldOfView> crops,
       pgmlink::ConsTrackingInferenceModel::IndexType numWeights,
       pgmlink::ConsTrackingInferenceModel::LabelType numLabels,
       int ndim,
       boost::shared_ptr<pgmlink::HypothesesGraph> hypothesesGraph
   )
   {
       std::cout << "numModels = " << numModels << " crops.size() = " << crops.size() << std::endl;
       if(numModels!=crops.size()){
           std::cout << "Number of crops and crops size do not match!!!" << std::endl;
           return;
       }

       using namespace pgmlink;

       this->lossParams_.resize(numModels);
       this->isCached_.resize(numModels);
       this->count_.resize(numModels,0);
       this->weights_ = Weights(numWeights);

       this->gms_.resize(numModels);
       this->gmsWithLoss_.resize(numModels);
   }

   void setGraphicalModel(size_t modelIndex, size_t numVariables){
       this->gms_[modelIndex].addVariable(numVariables);
   }
};





} // datasets
} // opengm







#endif // OPENGM_DATASET_H

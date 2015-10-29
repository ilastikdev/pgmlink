#ifndef OPENGM_DATASET_H
#define OPENGM_DATASET_H

#include <opengm/learning/dataset/editabledataset.hxx>

#include "pgmlink/inferencemodel/structuredlearningtrackinginferencemodel.h"
#include "pgmlink/tracking.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/graph.h"
#include "pgmlink/field_of_view.h"

namespace opengm {
namespace datasets{

template<class GM, class LOSS>
class StructuredLearningTrackingDataset : public EditableDataset<GM,LOSS,GM>{
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
       //pgmlink::StructuredLearningTrackingInferenceModel::IndexType numModels,
       std::vector<pgmlink::FieldOfView> crops,
       //pgmlink::StructuredLearningTrackingInferenceModel::IndexType numWeights,
       //pgmlink::StructuredLearningTrackingInferenceModel::LabelType numLabels,
       //int ndim,
       //boost::shared_ptr<pgmlink::HypothesesGraph> hypothesesGraph,
       opengm::learning::Weights<double>& trackingWeights
   )
   {
       pgmlink::StructuredLearningTrackingInferenceModel::IndexType numModels =crops.size();
//     std::cout << "Number of crops and crops size do not match!!!" << std::endl; // use assert
//           return;
//       }

       this->lossParams_.resize(numModels);
       this->isCached_.resize(numModels);
       this->count_.resize(numModels,0);
       this->getWeights() = trackingWeights;

       this->gms_.resize(numModels);
       this->gts_.resize(numModels);
       this->gmsWithLoss_.resize(numModels);
   }

   void setGraphicalModel(size_t m, GMType model){
       this->gms_[m] = model;
   }

   void build_model_with_loss(size_t m){
       this->buildModelWithLoss(m);
   }

   GMType graphicalModel(size_t m){
       return this->gms_[m];
   }

   void resizeGTS(size_t m){

       this->gts_[m].resize(this->gms_[m].numberOfVariables(), 0);
   }

   void setGTS(size_t modelIndex, size_t labelIndex, LabelType label){
       this->gts_[modelIndex][labelIndex] = label;
   }

};
} // datasets
} // opengm







#endif // OPENGM_DATASET_H

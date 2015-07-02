#ifndef OPENGM_DATASET_H
#define OPENGM_DATASET_H

#include <opengm/learning/dataset/dataset.hxx>
#include "pgmlink/inferencemodel/constrackinginferencemodel.h"
#include "pgmlink/hypotheses.h"
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
   {
       std::cout << "StructuredLearningTrackingDataset DEFAULT constructor 0 param" << std::endl;
   }

   StructuredLearningTrackingDataset(
       pgmlink::ConsTrackingInferenceModel::IndexType numModels=0,
       //std::vector<pgmlink::FieldOfView> crops,
       pgmlink::ConsTrackingInferenceModel::IndexType numWeights=0,
       pgmlink::ConsTrackingInferenceModel::LabelType numLabels=0//,
       //boost::shared_ptr<pgmlink::HypothesesGraph> hypothesesGraph
   )
   {};
};

} // datasets
} // opengm







#endif // OPENGM_DATASET_H

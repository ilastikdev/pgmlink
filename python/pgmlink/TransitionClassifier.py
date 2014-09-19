import numpy

#inspired by opObjectClassification
def make_feature_array(feats, selected):

    featsMatrix= []
    for featname in selected:
        featurea = feats[featname]
        for value in featurea:
            ft = numpy.reshape(numpy.asarray(value),(-1,))
            featsMatrix.append(ft)
        
    featsMatrix = numpy.concatenate(featsMatrix)
    return featsMatrix

class TransitionClassifier:
    def __init__(self,classifier,selected_features):
        self.classifier = classifier
        self.selected_features = selected_features
        
    def getDistance(self,Traxel1,Traxel2):
        """
        returns probability and variance of transition from Traxel1 to Traxel2
        based on transition classifier (gaussian process classifier)
        """
        
        Features1 = make_feature_array(Traxel1.features,self.selected_features)
        Features2 = make_feature_array(Traxel2.features,self.selected_features)
        
        X = numpy.reshape(numpy.concatenate((Features1,Features2)),(1,-1))
        
        classifier_prob = self.classifier.predict_probabilities(X)[0,0]
        classifier_var  = self.classifier.predict_probabilities(X)[0,1]
        return classifier_prob,classifier_var
    
    def test(self):
        print "Test"
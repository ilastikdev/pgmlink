import numpy as np

class TransitionClassifier:
    def __init__(self,classifier, selected_features):
        self.classifier = classifier
        if len(selected_features) != 1 or selected_features[0] != 'SquaredDistance<com>':
            raise NotImplementedError, 'other features are not supported yet'
        
    @staticmethod
    def getFeatures(traxel1, traxel2):
        # only squared distances for now
        return np.linalg.norm( np.array(traxel1.features['com']) - np.array(traxel2.features['com']) )

    def predict(self, traxel1, traxel2):
        """
        returns probability and variance of transition from Traxel1 to Traxel2
        based on transition classifier (gaussian process classifier)
        """
        x = self.getFeatures(traxel1, traxel2)
        
        return self.classifier.predit_probabilities(x, with_variance=True)

import sys

sys.path.append("../.")

import unittest as ut
import pgmlink
# from ilastik.applets.tracking.conservation.transitionClassifierTraining import trainTransition
import trainClassifiers
from ilastik.applets.tracking.conservation.TransitionClassifier import TransitionClassifier


def mk_traxel(x, y, z, id, t, fs):
    traxel = pgmlink.Traxel()
    traxel.ID = id
    traxel.Timestep = t
    traxel.set_feature_store(fs)

    traxel.add_feature_array("com", 3)
    traxel.set_feature_value("com", 0, x)
    traxel.set_feature_value("com", 1, y)
    traxel.set_feature_value("com", 2, z)
    
    traxel.add_feature_array("divProb",2)
    traxel.set_feature_value("divProb",0, 0.1)
    traxel.set_feature_value("divProb",1, 0.9)

    traxel.add_feature_array("divProb_Var", 1)
    traxel.set_feature_value("divProb_Var", 0, 0.01)

    traxel.add_feature_array("detProb", 2)
    traxel.set_feature_value("detProb", 0, 0.1)
    traxel.set_feature_value("detProb", 1, 0.9)

    traxel.add_feature_array("detProb_Var", 1)
    traxel.set_feature_value("detProb_Var", 0, 0.01)
    
    traxel.add_feature_array("count",1)
    traxel.set_feature_value("count",0, 1)
    return traxel

class Test_Traxels( ut.TestCase ):
    def runTest( self ):
        fs = pgmlink.FeatureStore()

        t = mk_traxel(34,23,12,77, 0, fs)
        self.assertEqual(t.ID, 77)
        # self.assertEqual(len(t.features), 1)
        self.assertEqual(t.get_feature_value("com", 0), 34)
        self.assertEqual(t.get_feature_value("com", 1), 23)
        self.assertEqual(t.get_feature_value("com", 2), 12)

class Test_TraxelStore( ut.TestCase ):
    def test_pickle( self ):
        import cPickle
        fs = pgmlink.FeatureStore()

        t1 = mk_traxel(1,2,3,33,0,fs)
        t2 = mk_traxel(5,6,7,44,0,fs)
        ts = pgmlink.TraxelStore()
        ts.add(fs, t1)
        ts.add(fs, t2)

        saved = cPickle.dumps(ts)
        loaded = cPickle.loads(saved)


class Test_HypothesesGraph( ut.TestCase ):
    def test_graph_interface( self ):
        # exercise the interface

        g = pgmlink.HypothesesGraph()
        n1 = g.addNode(0)
        n2 = g.addNode(5)
        n3 = g.addNode(7)
        a1 = g.addArc(n1, n2)
        a2 = g.addArc(n1,n3)
        self.assertEqual( pgmlink.countNodes(g), 3)
        self.assertEqual( pgmlink.countArcs(g), 2)
        self.assertEqual( g.earliest_timestep(), 0 )
        self.assertEqual( g.latest_timestep(), 7 )

        g.erase(a2)
        g.erase(n3)
        self.assertEqual( pgmlink.countNodes(g), 2)
        self.assertEqual( pgmlink.countArcs(g), 1)
        self.assertTrue( g.valid(n1) )
        self.assertTrue( g.valid(n2) )
        self.assertTrue( not g.valid(n3) )
        self.assertTrue( g.valid(a1) )
        self.assertTrue( not g.valid(a2) )

    def test_property_maps( self ):
        t = pgmlink.Traxel()
        t.Id = 33

        g = pgmlink.HypothesesGraph()
        n1 = g.addNode(0)
        m = g.addNodeTraxelMap()
        m[n1] = t
        self.assertEqual(m[n1].Id, 33)


class Test_CrossCorrelation( ut.TestCase ):
    def runTest( self ):
        pass
#         import numpy as np
#         img1 = np.array( [ [1, 0, 1], [1, 1, 1], [1, 1, 1] ], dtype=np.float)
#         img2 = np.array( [ [0, 0, 0], [1, 0, 1], [1, 1, 1] ], dtype=np.float)
# #        pgmlink.patchedCrossCorrelation(img1,img2,int(3),int(0),int(0),True)
#         print 'bla'
#         pgmlink.patchedCrossCorrelation(int(3),int(0),int(0),True)
#         print 'success'


def loadGPClassifier(fn, h5_group='/TransitionGPClassifier/'):
    import h5py
    from lazyflow.classifiers.gaussianProcessClassifier import GaussianProcessClassifier

    with h5py.File(fn, 'r') as f:
        g = f[h5_group]
        gpc = GaussianProcessClassifier()
        gpc.deserialize_hdf5(g)

        features = []
        for op in g['Features'].keys():
            for feat in g['Features'][op]:
                features.append('%s<%s>' % (op, feat))

    return gpc, features


class TestTransitionClassifier(ut.TestCase):

    def runTest( self ):
        print 'setting up traxels'
        ts = pgmlink.TraxelStore()
        fs = pgmlink.FeatureStore()

        t1 = mk_traxel(1,2,3,33,0,fs)
        t2 = mk_traxel(5,6,7,44,1,fs)

        ts.add(fs, t1)
        ts.add(fs, t2)

        sigmas = pgmlink.VectorOfDouble()
        sigmas.append(0) # appearance variance
        sigmas.append(0) # disappearance variance

        distr_id = pgmlink.DistrId.ClassifierUncertainty
        uncertainty_parameter = pgmlink.UncertaintyParameter(3, # number of iterations
                                                             distr_id, # distribution id
                                                             sigmas)

        fov = pgmlink.FieldOfView(0,
                                  0,
                                  0,
                                  0,
                                  2,
                                  2,
                                  2,
                                  2,)

        print 'initialize conservation tracking'
        tracker = pgmlink.ConsTracking(
                                1, # max_num_objects
                                True, # size_dependent_detection_prob
                                1, # average object size
                                10,  # max neighbor distance
                                False, # with divisions
                                0, # division threshold
                                "none", # rf filename
                                fov, # field of view
                                "none", # dump traxelstore
                                )

        print 'build hypotheses graph'
        tracker.buildGraph(ts)

        print 'load pre-trained GP classifier'
        gpc, selected_features = loadGPClassifier(fn='./test-transition-classifier.h5', h5_group='/TransitionGPClassifier/')
        classifier = TransitionClassifier(gpc, selected_features)

        print 'track...'
        tracker.track(
            0., # forbidden costs
            0.01, # ep_gap
            False, # with tracklets
            10., # detection weight
            10., # division weight
            10., # transition weight
            100., # disappearance cost
            100., # appearance cost
            False, # with merger resolution
            3, # ndim
            10., # transition parameter
            0., # border width
            True, # with constraints
            uncertainty_parameter,
            1e+75, # cplex timeout
            classifier, # transition classifier
        )

        print 'done.'
        
if __name__=="__main__":
    ut.main()


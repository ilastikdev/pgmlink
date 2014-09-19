import sys
sys.path.append("../.")

import unittest as ut
import pgmlink

def mk_traxel(x,y,z,id):
    t = pgmlink.Traxel()
    t.ID = id
    t.add_feature_array("com", 3)
    t.set_feature_value("com", 0, x)
    t.set_feature_value("com", 1, y)
    t.set_feature_value("com", 2, z)
    
    t.add_feature_array("divProb",2)
    t.set_feature_value("divProb",0, 0.1)
    t.set_feature_value("divProb",1, 0.9)
    
    t.add_feature_array("count",1)
    t.set_feature_value("count",0, 1)
    return t

class Test_Traxels( ut.TestCase ):
    def runTest( self ):
        t = mk_traxel(34,23,12,77)
        self.assertEqual(t.ID, 77)
        self.assertEqual(len(t.features), 1)
        self.assertEqual(t.get_feature_value("com", 0), 34)
        self.assertEqual(t.get_feature_value("com", 1), 23)
        self.assertEqual(t.get_feature_value("com", 2), 12)

class Test_TraxelStore( ut.TestCase ):
    def test_pickle( self ):
        import cPickle
        t1 = mk_traxel(1,2,3,33)
        t2 = mk_traxel(5,6,7,44)
        ts = pgmlink.TraxelStore()
        ts.add(t1)
        ts.add(t2)

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
        import numpy as np
        img1 = np.array( [ [1, 0, 1], [1, 1, 1], [1, 1, 1] ], dtype=np.float)
        img2 = np.array( [ [0, 0, 0], [1, 0, 1], [1, 1, 1] ], dtype=np.float)
#        pgmlink.patchedCrossCorrelation(img1,img2,int(3),int(0),int(0),True)
        print 'bla'
        pgmlink.patchedCrossCorrelation(int(3),int(0),int(0),True)
        print 'success'
class TestTransitionClassifier(ut.TestCase):
    def runTest( self ):
        ts = pgmlink.TraxelStore()
        t1 = mk_traxel(1,2,3,33)
        t2 = mk_traxel(5,6,7,44)
        ts = pgmlink.TraxelStore()
        ts.add(t1)
        ts.add(t2)
        
        distr = pgmlink.DistrId.GaussianPertubation
        
        up = pgmlink.UncertaintyParameter(1,distr,pgmlink.VectorOfDouble())
        fov = pgmlink.FieldOfView(0,
                                  0,
                                  0,
                                  0,
                                  2,
                                  2,
                                  2,
                                  2,)
         
        tracker = pgmlink.ConsTracking(1,
                                1,
                                1,
                                "none",  # detection_rf_filename
                                1,   # size_dependent_detection_prob
                                0,       # forbidden_cost
                                1., # ep_gap
                                1., # median_object_size
                                False, #withTracklets,
                                1.,
                                1.,
                                True,
                                1., # disappearance cost
                                1., # appearance cost
                                False,
                                4,
                                1.,
                                1.,
                                fov,
                                True, #with_constraints
                                up,
                                1.,
                                "none" # dump traxelstore
                                )
        cm = pgmlink.TimestepIdCoordinateMap()
        cm.initialize()
        classifier = pgmlink.TransitionClassifier(None,["con"])
        tracker(ts,cm.get(),classifier)
        
if __name__=="__main__":
    ut.main()

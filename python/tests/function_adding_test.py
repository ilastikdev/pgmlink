import sys

sys.path.append("../.")

import unittest as ut
import pgmlink
import numpy as np

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

def my_detection_func(traxel, state):
    prob = traxel.get_feature_value("detProb", state)
    return 10.0*(1.0-prob)

def my_transition_func(distance):
    return 10.0*float(distance)**2

def my_appearance_func(traxel):
    pos = np.array([traxel.get_feature_value("com", i) for i in range(3)])
    return np.linalg.norm(pos)

class TestFunctionAdding(ut.TestCase):

    def runTest( self ):
        print 'setting up traxels'
        ts = pgmlink.TraxelStore()
        fs = pgmlink.FeatureStore()

        t1 = mk_traxel(1,2,3,33,0,fs)
        t2 = mk_traxel(5,6,7,44,1,fs)

        ts.add(fs, t1)
        ts.add(fs, t2)

        sigmas = pgmlink.VectorOfDouble()
        uncertaintyParams = pgmlink.UncertaintyParameter(1, pgmlink.DistrId.PerturbAndMAP, sigmas)

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
                                "none", # dump event vector
                                pgmlink.ConsTrackingSolverType.CplexSolver
                                )

        print 'build hypotheses graph'
        tracker.buildGraph(ts)

        print 'setup parameters'
        params = tracker.get_conservation_tracking_parameters(
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
            uncertaintyParams,
            1e+75, # cplex timeout
            None, # TransitionClassifier
            pgmlink.ConsTrackingSolverType.CplexSolver, # Solver
            1,  # number of threads for CPLEX
            False) # whether or not to use an additional cross-timestep constraint

        print 'set special functions'
        params.register_detection_func(my_detection_func)
        params.register_transition_func(my_transition_func)
        params.register_appearance_func(my_appearance_func)

        print 'track...'
        tracker.track(params, False)

        print 'done.'
        
if __name__=="__main__":
    ut.main()


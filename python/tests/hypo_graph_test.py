import pgmlink
import math
from random import choice

import sys
sys.path.append("../.")

import unittest as ut
import pgmlink

def mk_traxel(x,y,z,tid,t):
    tr = pgmlink.Traxel() 
    tr.Id = tid
    tr.Timestep = t
    tr.add_feature_array("com", 3)
    tr.set_feature_value("com", 0, x)
    tr.set_feature_value("com", 1, y)
    tr.set_feature_value("com", 2, z)
    return tr

        
if __name__=="__main__":
    
    t1 = mk_traxel(1,2,3,33,1)
    t2 = mk_traxel(5,6,7,44,2)
    t3 = mk_traxel(6,6,6,55,3)
    ts = pgmlink.TraxelStore()
    ts.add(t1)
    ts.add(t2)
    ts.add(t3)

    ###########################################
    
    print "initialize tracking and obtain graph"
    
    fov = pgmlink.FieldOfView(0,
                                      0 ,
                                      0 ,
                                      0 ,
                                      3 ,
                                      1000 ,
                                      1000 ,
                                      1000 )

    tracker = pgmlink.ConsTracking(3,False,30.,20,False,0.3,"none",fov,"none")

    g = tracker.buildGraph(ts)
    
    print "#nodes\t",pgmlink.countNodes(g)
    print "#arcs\t",pgmlink.countArcs(g)
    
    ###########################################
    
    print "add Nodes with Traxel"
        
    new_node_list = []
    
    new_node_list.append(g.addTraxel(mk_traxel(2,2,2,66,1)))
    new_node_list.append(g.addTraxel(mk_traxel(3,3,3,77,2)))
    new_node_list.append(g.addTraxel(mk_traxel(4,4,4,88,3)))

    print "#nodes\t",pgmlink.countNodes(g)
    print "#arcs\t",pgmlink.countArcs(g)

    ###########################################
        
    print "connect all new nodes"
    
    g.addArc(new_node_list[0],new_node_list[1])
    g.addArc(new_node_list[1],new_node_list[2])
    
    print "#nodes\t",pgmlink.countNodes(g)
    print "#arcs\t",pgmlink.countArcs(g)

    ###########################################
    
    print "label graph randomly and build dictionary"
    
    n_it = pgmlink.NodeIt(g)
    a_it = pgmlink.ArcIt(g)

    nodeLabels = {}
    arcLabels = {}
    
    tmap = g.getNodeTraxelMap()
    
    for i in range(g.earliest_timestep(), g.latest_timestep()+1):
        nodeLabels[i] = {}
        arcLabels[i] = {}
    
    for n in n_it:
        random_label = choice([True,False])
        nodeLabels[tmap[n].Timestep][tmap[n].Id] = random_label
        g.addNodeLabel(n,random_label) 
    
    for a in a_it:
        random_label = choice([True,False])
        arcLabels[tmap[g.source(a)].Timestep][tmap[g.source(a)].Id] = dict([(tmap[g.target(a)].Id, choice([True,False]))])
        g.addArcLabel(a,random_label)
        
    print nodeLabels
    print arcLabels       
    

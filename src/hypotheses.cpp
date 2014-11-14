#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <algorithm>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/tuple/tuple.hpp>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/nearest_neighbors.h"
#include "pgmlink/traxels.h"

using namespace std;

namespace pgmlink {
////
//// class HypothesesGraph
////
HypothesesGraph::Node HypothesesGraph::add_node(node_timestep_map::Value timestep) {
    node_timestep_map& timestep_m = get(node_timestep());

    HypothesesGraph::Node node = addNode();
    timestep_m.set(node, timestep);
    timesteps_.insert( timestep );
    return node;
}

HypothesesGraph::Node HypothesesGraph::add_node(std::vector<node_timestep_map::Value> timesteps) {
    node_timestep_map& timestep_m = get(node_timestep());

    HypothesesGraph::Node node = addNode();
    for (std::vector<node_timestep_map::Value>::const_iterator it = timesteps.begin(); it!=timesteps.end(); ++it) {
        timestep_m.set(node, *it);
        timesteps_.insert( *it );
    }
    return node;
}

void HypothesesGraph::init_labeling_maps() {
    add(appearance_label()).add(disappearance_label()).add(division_label()).add(arc_label());

    property_map< appearance_label, HypothesesGraph::base_graph>::type& map1 = get(appearance_label());
    property_map< disappearance_label, HypothesesGraph::base_graph>::type& map2 = get(disappearance_label());
    property_map< division_label, HypothesesGraph::base_graph>::type& map3 = get(division_label());
    for (HypothesesGraph::NodeIt it(*this); it != lemon::INVALID; ++it) {
        map1.set(it, 0);
        map2.set(it, 0);
        map3.set(it, 0);
    }

    property_map< arc_label, HypothesesGraph::base_graph>::type& map = get(arc_label());
    for (HypothesesGraph::ArcIt it(*this); it != lemon::INVALID; ++it) {
        map.set(it, 0);
    }
}

void HypothesesGraph::add_appearance_label(HypothesesGraph::Node node,label_type label) {  
  property_map< appearance_label, HypothesesGraph::base_graph>::type& gt_label = get(appearance_label());
  gt_label.set(node,label);
}

void HypothesesGraph::add_disappearance_label(HypothesesGraph::Node node,label_type label) {  
  property_map< disappearance_label, HypothesesGraph::base_graph>::type& gt_label = get(disappearance_label());
  gt_label.set(node,label);
}

void HypothesesGraph::add_division_label(HypothesesGraph::Node node,label_type label) {
  property_map< division_label, HypothesesGraph::base_graph>::type& gt_label = get(division_label());
  gt_label.set(node,label);
}

void HypothesesGraph::add_arc_label(HypothesesGraph::Arc arc,label_type label) {
  property_map<  arc_label, HypothesesGraph::base_graph>::type& gt_label = get(arc_label());
  gt_label.set(arc,label);
}

HypothesesGraph::Node HypothesesGraph::add_traxel(Traxel ts){
  LOG(logDEBUG4) << "add traxel(id="<<ts.Id<<") at t="<< ts.Timestep;
  HypothesesGraph::Node node = add_node(ts.Timestep);
  get(node_traxel()).set(node,ts);
  return node;
}

  const std::set<HypothesesGraph::node_timestep_map::Value>& HypothesesGraph::timesteps() const {
    return timesteps_;
}

HypothesesGraph::node_timestep_map::Value HypothesesGraph::earliest_timestep() const {
    return *(timesteps_.begin());
}

HypothesesGraph::node_timestep_map::Value HypothesesGraph::latest_timestep() const {
    return *(timesteps_.rbegin());
}



HypothesesGraph& prune_inactive(HypothesesGraph& g) {
	  LOG(logDEBUG) << "prune_inactive(): entered";
	  property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());

	  property_map<node_active, HypothesesGraph::base_graph>::type* active_nodes = 0;
    property_map<node_active2, HypothesesGraph::base_graph>::type* active2_nodes = 0;
    bool active2_used = false;
    if (g.getProperties().count("node_active") > 0) {
		  active_nodes = &g.get(node_active());
    } else {
		  assert(g.getProperties().count("node_active2") > 0);
		  active2_nodes = &g.get(node_active2());
		  active2_used = true;
    }

	  // prune inactive arcs
	  LOG(logDEBUG) << "prune_inactive: prune inactive arcs";
	  typedef property_map<arc_active, HypothesesGraph::base_graph>::type::FalseIt inactive_arc_it;

	  // we first collect and then erase (iterator will be invalid after erase; therefore a two-step
	  // procedure)

	  // collect inactive arcs
	  vector<HypothesesGraph::Arc> arcs_to_prune;

	  for(inactive_arc_it it(active_arcs); it!=lemon::INVALID; ++it) {
		  arcs_to_prune.push_back(it);
		  assert(g.valid(it));
		  LOG(logDEBUG1) << "prune_inactive: arc to be pruned: " << g.id(it);
	  }

	  std::sort(arcs_to_prune.begin(), arcs_to_prune.end());
	  std::reverse(arcs_to_prune.begin(), arcs_to_prune.end());

	  // prune inactive arcs
	  for(vector<HypothesesGraph::Arc>::const_iterator it = arcs_to_prune.begin(); it!= arcs_to_prune.end(); ++it) {
		  if (g.valid(*it)) {
			  LOG(logDEBUG1) << "prune_inactive: pruned arc: " << g.id(*it);
			  g.erase(*it);
		  }
		  assert(!g.valid(*it));
	  }

	  // prune inactive nodes
	  LOG(logDEBUG) << "prune_inactive: prune inactive nodes";
	  // collect inactive nodes
	  vector<HypothesesGraph::Node> nodes_to_prune;

	  if (active2_used) {
		  for (HypothesesGraph::NodeIt it(g); it != lemon::INVALID; ++it) {
			  if (!(*active2_nodes)[it]) {
				  nodes_to_prune.push_back(it);
				  assert(g.valid(it));
			  }
		  }
	  } else {
		  typedef property_map<node_active, HypothesesGraph::base_graph>::type::FalseIt inactive_node_it;
		  for (inactive_node_it it(*active_nodes); it != lemon::INVALID; ++it) {
			  nodes_to_prune.push_back(it);
			  assert(g.valid(it));
		  }
	  }

	  property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
	  // prune inactive nodes
	  for(vector<HypothesesGraph::Node>::const_iterator it = nodes_to_prune.begin(); it!= nodes_to_prune.end(); ++it) {
		  LOG(logDEBUG3) << "prune_inactive: prune node: " << g.id(*it) << ", Traxel = " << traxel_map[*it];
		  for(HypothesesGraph::OutArcIt arcit(g,*it); arcit != lemon::INVALID; ++arcit) {
			  assert(active_arcs[arcit] == false && "arc may not be active");
		  }
		  for(HypothesesGraph::InArcIt arcit(g,*it); arcit != lemon::INVALID; ++arcit) {
			  assert(active_arcs[arcit] == false && "arc may not be active");
		  }
		  g.erase(*it);
		  assert(!g.valid(*it));
	  }
	  LOG (logDEBUG) << "prune_inactive(): done";
	  return g;
}


// Introducing: Three awkward workarounds.
// We don't want to have ifs everywhere,
// so we decided to write these three magnificent convenience functions
size_t get_active_node(
		property_map<node_active, HypothesesGraph::base_graph>::type* nodes,
		property_map<node_active2, HypothesesGraph::base_graph>::type* nodes2,
		property_map<node_active_count, HypothesesGraph::base_graph>::type* nodes_vector,
		size_t state,
		HypothesesGraph::Node node,
		size_t iterationStep=0){
	if (state ==0){
		return (*nodes)[node];
	}
	else if (state==1){
		return (*nodes2)[node];
	}
	return (*nodes_vector)[node][iterationStep];
}

bool get_active_arc(
		property_map<arc_active, HypothesesGraph::base_graph>::type* arcs,
		property_map<arc_active_count, HypothesesGraph::base_graph>::type* arcs_vector,
		size_t state,
		HypothesesGraph::Arc arc,
		size_t iterationStep=0){
	if (state==2){
		return (*arcs_vector)[arc][iterationStep];
	}
	return (*arcs)[arc];
}

bool get_active_division(
		property_map<division_active, HypothesesGraph::base_graph>::type* divs,
		property_map<division_active_count, HypothesesGraph::base_graph>::type* divs_vector,
		size_t state,
		HypothesesGraph::Node node,
		size_t iterationStep=0){
	if (state==2){
		return (*divs_vector)[node][iterationStep];
	}
	return (*divs)[node];
}

boost::shared_ptr<std::vector< std::vector<Event> > > events(const HypothesesGraph& g, int iterationStep) {

	LOG(logDEBUG) << "events(): entered";
	boost::shared_ptr<std::vector< std::vector<Event> > > ret(new vector< vector<Event> >);
	typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
	node_timestep_map_t& node_timestep_map = g.get(node_timestep());
	typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map_t;
	node_traxel_map_t& node_traxel_map = g.get(node_traxel());

	property_map<node_active, HypothesesGraph::base_graph>::type* nodes=0;
	property_map<node_active2, HypothesesGraph::base_graph>::type* nodes2=0;
	property_map<node_active_count, HypothesesGraph::base_graph>::type* nodes_vector=0;


	property_map<arc_active, HypothesesGraph::base_graph>::type* arcs=0;
	property_map<arc_active_count, HypothesesGraph::base_graph>::type* arcs_vector=0;

	bool with_division_detection = false;
	property_map<division_active, HypothesesGraph::base_graph>::type* divisions=0;
	property_map<division_active_count, HypothesesGraph::base_graph>::type* divisions_vector=0;

	size_t vector_state;
	bool with_mergers = false;
	if (g.getProperties().count("node_active_count") > 0){
		vector_state = 2;
		arcs_vector=  &g.get(arc_active_count());
		nodes_vector = &g.get(node_active_count());
		if (g.getProperties().count("node_active2") > 0) {
			with_mergers = true;}
	}
	else if (g.getProperties().count("node_active2") > 0) {
		with_mergers = true;
		nodes2 = &g.get(node_active2());
		arcs=  &g.get(arc_active());
		LOG(logDEBUG1) << "events(): with_mergers = true";
		vector_state = 1;
	}
	else {
		vector_state = 0;
		arcs=  &g.get(arc_active());
		nodes = &g.get(node_active());
	}

	if (g.getProperties().count("division_active") > 0) {
		with_division_detection = true;
		if (vector_state==2) {
			divisions_vector = &g.get(division_active_count());
		} else {
			divisions = &g.get(division_active());
		}
	}

	bool with_origin = false;
	property_map<node_originated_from, HypothesesGraph::base_graph>::type* origin_map;
	if (g.getProperties().count("node_originated_from") > 0) {
		origin_map = &g.get(node_originated_from());
		with_origin = true;
		LOG(logDEBUG1) << "events(): with_origin enabled";
	}

	// for every timestep
	LOG(logDEBUG1) << "events(): earliest_timestep: " << g.earliest_timestep();
	LOG(logDEBUG1) << "events(): latest_timestep: " << g.latest_timestep();

	// add an empty first timestep
	ret->push_back(vector<Event>());

	for(int t = g.earliest_timestep(); t < g.latest_timestep(); ++t) {
		LOG(logDEBUG2) << "events(): processing timestep: " << t;
		ret->push_back(vector<Event>());

		map<unsigned int, vector<unsigned int> > resolver_map;

		// for every node: destiny
		LOG(logDEBUG2) << "events(): for every node: destiny";
		for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t); node_at!=lemon::INVALID; ++node_at) {
			assert(node_traxel_map[node_at].Timestep == t);
			LOG(logDEBUG4) <<t<< " "<< node_traxel_map[node_at].Id;

			if(!get_active_node(nodes,nodes2,nodes_vector,vector_state,node_at,iterationStep)){
				continue;
			}

			if (with_origin && (*origin_map)[node_at].size() > 0 && t > g.earliest_timestep()) {
				LOG(logINFO) << "events(): collecting resolver node ids for all merger nodes " << t << ", " << (*origin_map)[node_at][0];
				resolver_map[(*origin_map)[node_at][0]].push_back(node_traxel_map[node_at].Id);
				const std::vector<float>& tmp_feat = (node_traxel_map[node_at].features.find("com"))->second; //node_traxel_map[node_at].features["com"];
				std::copy(tmp_feat.begin(), tmp_feat.end(),
						std::back_insert_iterator<std::vector<unsigned> >(resolver_map[(*origin_map)[node_at][0]]));
			}

			LOG(logDEBUG3) << "Number of detected objects: " << get_active_node(nodes,nodes2,nodes_vector,vector_state,node_at,iterationStep);

			// count outgoing arcs

			size_t count = 0;
			for(HypothesesGraph::base_graph::OutArcIt a(g, node_at); a!=lemon::INVALID; ++a) {
				if (get_active_arc(arcs,arcs_vector,vector_state,a,iterationStep)){
					++count;
				}
			}

			LOG(logDEBUG3) << "events(): counted outgoing arcs: " << count;

			// construct suitable Event object
			switch(count) {
			// Disappearance
			case 0: {
				if (t<g.latest_timestep()) {
					Event e;
					e.type = Event::Disappearance;
					e.traxel_ids.push_back(node_traxel_map[node_at].Id);
					(*ret)[t-g.earliest_timestep()+1].push_back(e);
					LOG(logDEBUG3) << e;
				}
				break;
			}
			// Move
			case 1: {
				Event e;
				e.type = Event::Move;
				e.traxel_ids.push_back(node_traxel_map[node_at].Id);

				for(HypothesesGraph::base_graph::OutArcIt a(g, node_at); a != lemon::INVALID; ++a) {
					if (get_active_arc(arcs,arcs_vector,vector_state,a,iterationStep)){
						e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
					}

				}
				(*ret)[t-g.earliest_timestep()+1].push_back(e);

				LOG(logDEBUG3) << e;
				break;
			}
			// Division or Splitting
			default: {
				Event e;
				if (with_division_detection) {
					if (count == 2 && get_active_division(divisions,divisions_vector,vector_state,node_at,iterationStep)) {
						e.type = Event::Division;
						e.traxel_ids.push_back(node_traxel_map[node_at].Id);
						for(HypothesesGraph::base_graph::OutArcIt a(g, node_at); a != lemon::INVALID; ++a) {
							if (get_active_arc(arcs,arcs_vector,vector_state,a,iterationStep)){
								e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
							}
						}
						(*ret)[t-g.earliest_timestep()+1].push_back(e);


						LOG(logDEBUG3) << e;
					} else {

						for(HypothesesGraph::base_graph::OutArcIt a(g, node_at); a != lemon::INVALID; ++a) {
							if (get_active_arc(arcs,arcs_vector,vector_state,a,iterationStep)) {

								e.type = Event::Move;
								e.traxel_ids.clear();
								e.traxel_ids.push_back(node_traxel_map[node_at].Id);
								e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
								(*ret)[t-g.earliest_timestep()+1].push_back(e);
								LOG(logDEBUG3) << e;
							}
						}
					}
				} else { // for backward compatibility
					if (count != 2) {
						throw runtime_error("events(): encountered node dividing in three or more nodes in graph");
					}
					e.type = Event::Division;
					e.traxel_ids.push_back(node_traxel_map[node_at].Id);
					for(HypothesesGraph::base_graph::OutArcIt a(g, node_at); a != lemon::INVALID; ++a) {
						if (get_active_arc(arcs,arcs_vector,vector_state,a,iterationStep)){
							e.traxel_ids.push_back(node_traxel_map[g.target(a)].Id);
						}
					}
					(*ret)[t-g.earliest_timestep()+1].push_back(e);

					LOG(logDEBUG3) << e;
				}
				break;
			}
			}
		}
		for (map<unsigned int, vector<unsigned int> >::iterator map_it = resolver_map.begin(); map_it != resolver_map.end(); ++map_it) {
			Event e;
			e.type = Event::ResolvedTo;
			e.traxel_ids.push_back(map_it->first);
			for (std::vector<unsigned int>::iterator it = map_it->second.begin(); it != map_it->second.end(); ++it) {
				e.traxel_ids.push_back(*it);
			}
			(*ret)[t-g.earliest_timestep()-1].push_back(e);
			LOG(logDEBUG1) << e;
		}


		// appearances and mergers in next timestep
		LOG(logDEBUG2) << "events(): appearances in next timestep";
		if (t+1 <= g.latest_timestep()) {
			for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t+1); node_at!=lemon::INVALID; ++node_at) {
				// count incoming arcs
				if (!get_active_node(nodes,nodes2,nodes_vector,vector_state,node_at,iterationStep)){
					continue;
				}
				int count = 0;
				for(HypothesesGraph::base_graph::InArcIt a(g, node_at); a!=lemon::INVALID; ++a) {
					if (get_active_arc(arcs,arcs_vector,vector_state,a,iterationStep)) ++count;
				}
				LOG(logDEBUG3) << "events(): counted incoming arcs in next timestep: " << count;
				// no incoming arcs => appearance
				if(count == 0 && t + 1 > g.earliest_timestep()) {
					Event e;
					e.type = Event::Appearance;
					e.traxel_ids.push_back(node_traxel_map[node_at].Id);
					(*ret)[t-g.earliest_timestep()+1].push_back(e);
					LOG(logDEBUG3) << e;
				}
			}
		}
		for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t); node_at!=lemon::INVALID; ++node_at) {
			if(with_mergers && get_active_node(nodes,nodes2,nodes_vector,vector_state,node_at,iterationStep) > 1) {
				Event e;
				e.type = Event::Merger;
				e.traxel_ids.push_back(node_traxel_map[node_at].Id);
				e.traxel_ids.push_back(get_active_node(nodes,nodes2,nodes_vector,vector_state,node_at,iterationStep));
				(*ret)[t-g.earliest_timestep()].push_back(e);
				LOG(logDEBUG3) << e;
			}
		}
	}

	LOG(logDEBUG2) << "events(): last timestep: " << g.latest_timestep();
	for(node_timestep_map_t::ItemIt node_at(node_timestep_map, g.latest_timestep()); node_at!=lemon::INVALID; ++node_at) {
		if(with_mergers && get_active_node(nodes,nodes2,nodes_vector,vector_state,node_at,iterationStep)> 1) {
			Event e;
			e.type = Event::Merger;
			e.traxel_ids.push_back(node_traxel_map[node_at].Id);
			e.traxel_ids.push_back(get_active_node(nodes,nodes2,nodes_vector,vector_state,node_at,iterationStep));
			(*ret)[g.latest_timestep()-g.earliest_timestep()].push_back(e);
			LOG(logDEBUG3) << e;
		}
	}
	LOG(logDEBUG2) << "events(): done.";
	return ret;
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& values)
{
    stream << "(";
    for(T t : values)
        stream << t << ", ";
    stream << ")";
    return stream;
}

boost::shared_ptr<std::vector< std::vector<Event> > > multi_frame_move_events(const HypothesesGraph& g) {

    boost::shared_ptr<std::vector< std::vector<Event> > > ret(new vector< vector<Event> >);
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    node_timestep_map_t& node_timestep_map = g.get(node_timestep());
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map_t;
    node_traxel_map_t& node_traxel_map = g.get(node_traxel());
    typedef property_map<node_originated_from, HypothesesGraph::base_graph>::type origin_map_t;
    origin_map_t& origin_map = g.get(node_originated_from());

    std::map<int, std::vector<Event> > multi_frame_move_map;

    // add an empty first timestep
    ret->push_back(vector<Event>());

    for(int t = g.earliest_timestep(); t < g.latest_timestep(); ++t) {
        LOG(logDEBUG2) << "events(): processing timestep: " << t;
        ret->push_back(vector<Event>());
        for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t); node_at!=lemon::INVALID; ++node_at) {
            assert(node_traxel_map[node_at].Timestep == t);
            if (origin_map[node_at].size()) {
                for(HypothesesGraph::InArcIt in_it(g, node_at); in_it != lemon::INVALID; ++in_it) {
                    HypothesesGraph::Node src_node = g.source(in_it);
                    if (origin_map[src_node].size()) {
                        break;
                    }
                    Event e;
                    e.type = Event::MultiFrameMove;
                    Traxel trax = node_traxel_map[src_node];
                    e.traxel_ids.push_back(trax.Id);
                    int t_local = t+1;
                    HypothesesGraph::Node n = node_at;
                    while (t_local <= g.latest_timestep()) {
                        HypothesesGraph::OutArcIt merge_it(g, n);
                        if (merge_it == lemon::INVALID) {
                            break;
                        }
                        n = g.target(merge_it);
                        assert(t_local == node_timestep_map[n]);
                        if (!origin_map[n].size()) {
                            trax = node_traxel_map[n];
                            assert(t_local == trax.Timestep);
                            e.traxel_ids.push_back(trax.Id);
                            e.traxel_ids.push_back(t-1-g.earliest_timestep());
                            multi_frame_move_map[t_local-g.earliest_timestep()].push_back(e);
                            break;
                        }
                        ++t_local;
                    }
                }
            }
        }

    }

    std::map<int, std::vector<Event> >::iterator map_it = multi_frame_move_map.begin();
    for (; map_it != multi_frame_move_map.end(); ++map_it) {
        (*ret)[map_it->first] = map_it->second;
    }
    
    return ret;
} /* multi_frame_move_events */

boost::shared_ptr< EventVectorVector > resolved_to_events(const HypothesesGraph& g)
{
    boost::shared_ptr<std::vector< std::vector<Event> > > ret(new vector< vector<Event> >);

    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    node_timestep_map_t& node_timestep_map = g.get(node_timestep());
    typedef property_map<merger_resolved_to, HypothesesGraph::base_graph>::type merger_resolved_map_t;
    merger_resolved_map_t& merger_resolved_map = g.get(merger_resolved_to());
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map_t;
    node_traxel_map_t& node_traxel_map = g.get(node_traxel());

    for(int t = g.earliest_timestep(); t <= g.latest_timestep(); ++t) {
        LOG(logDEBUG2) << "events(): processing timestep: " << t;
        ret->push_back(vector<Event>());
        for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t); node_at!=lemon::INVALID; ++node_at) {
            if(merger_resolved_map[node_at].size() > 0)
            {
                std::vector<unsigned int> new_traxel_ids = merger_resolved_map[node_at];
                Event e;
                e.type = Event::ResolvedTo;

                // insert original node
                Traxel trax = node_traxel_map[node_at];
                e.traxel_ids.push_back(trax.Id);

                // insert new nodes
                e.traxel_ids.insert(e.traxel_ids.end(), new_traxel_ids.begin(), new_traxel_ids.end());

                // insert into returned vector
                ret->back().push_back(e);

                LOG(logDEBUG2) << "Inserting Event::ResolvedTo at timestep " << t << " for traxel: " << trax.Id << " into: " << new_traxel_ids;
            }
        }
    }

    return ret;
}

boost::shared_ptr<std::vector< std::vector<Event> > > merge_event_vectors(const std::vector<std::vector<Event> >& ev1, const std::vector<std::vector<Event> >& ev2) {

	assert(ev1.size() == ev2.size());
    boost::shared_ptr<std::vector< std::vector<Event> > > ret(new vector< vector<Event> >);
    std::vector<std::vector<Event> >::const_iterator it1 = ev1.begin();
    std::vector<std::vector<Event> >::const_iterator it2 = ev2.begin();
    for (; it1 != ev1.end(); ++it1, ++it2) {
        ret->push_back(vector<Event>());
        std::back_insert_iterator<vector<Event> > push_back_inserter(*(ret->rbegin()));
        std::copy(it1->begin(), it1->end(), push_back_inserter);
        std::copy(it2->begin(), it2->end(), push_back_inserter);
    }
    return ret;
}

//
// generateTrackletGraph
//
// a traxel graph (containing nodes of traxels) with some active arcs is converted into a tracklet graph
// where all nodes connected by active paths are summarized in a single tracklet node
void generateTrackletGraph(const HypothesesGraph& traxel_graph, HypothesesGraph& tracklet_graph) {

    // go through the traxels graph, add each node which doesn't have an active incoming arc, and
    // follow the active outgoing path to add those nodes to the tracklet
    property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = traxel_graph.get(arc_active());
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = traxel_graph.get(node_traxel());

    property_map<node_tracklet, HypothesesGraph::base_graph>::type* traxels_tracklet_map;
    bool traxel_nodes_are_tracklets = false;
    if (traxel_graph.getProperties().count("node_tracklet") > 0) {
        traxel_nodes_are_tracklets = true;
        traxels_tracklet_map = &traxel_graph.get(node_tracklet());
    }

    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    node_timestep_map_t& node_timestep_map = traxel_graph.get(node_timestep());

    // add empty traxel_map to the tracklet graph in order to make the tracklet graph equivalent to traxelgraphs
    tracklet_graph.add(node_traxel());

    tracklet_graph.add(node_tracklet());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = tracklet_graph.get(node_tracklet());

    // this map stores the traxel_graph nodes at t+1 as keys which will be linked by the
    // list of tracklet nodes at t stored as map-values
    std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > node_link_map;

    // add nodes
    for(int t = traxel_graph.earliest_timestep(); t <= traxel_graph.latest_timestep(); ++t) {
        for(node_timestep_map_t::ItemIt traxel_node(node_timestep_map, t); traxel_node!=lemon::INVALID; ++traxel_node) {
            bool has_active_incoming = false;
            std::vector<HypothesesGraph::Arc> tracklet_incoming_arcs;
            for(HypothesesGraph::InArcIt a(traxel_graph, traxel_node); a != lemon::INVALID; ++a) {
                tracklet_incoming_arcs.push_back(a);
                if (active_arcs[a]) {
                    has_active_incoming = true;
                    break;
                }
            }

            if (has_active_incoming) {
                // if the traxel node has an active incoming arc, it has already been added to some tracklet
                continue;
            }
            std::vector<int> timesteps;
            std::vector<Traxel> tracklet;
            std::vector<Traxel> traxels;
            if (!traxel_nodes_are_tracklets) {
                traxels.push_back(traxel_map[traxel_node]);
            } else {
                traxels.insert(traxels.begin(), (*traxels_tracklet_map)[traxel_node].begin(),
                               (*traxels_tracklet_map)[traxel_node].end());
            }

            for (std::vector<Traxel>::const_iterator tr = traxels.begin(); tr != traxels.end(); ++tr) {
                timesteps.push_back(tr->Timestep);
                tracklet.push_back(*tr);
            }
            HypothesesGraph::Node tn = traxel_node;
            HypothesesGraph::Node tn_next;
            std::vector<HypothesesGraph::Arc> tn_outarcs;

            // follow the active outgoing path to add those nodes to the tracklet
            while (true) {
                bool active_outgoing = false;
                tn_outarcs.clear();

                for(HypothesesGraph::OutArcIt a(traxel_graph, tn); a != lemon::INVALID; ++a) {
                    tn_outarcs.push_back(a);

                    if (active_arcs[a]) {
                        assert(!active_outgoing); // "found more than one active outgoing arc"
                        tn_next = traxel_graph.target(a);

                        traxels.clear();
                        if (!traxel_nodes_are_tracklets) {
                            traxels.push_back(traxel_map[tn_next]);
                        } else {
                            traxels.insert(traxels.begin(), (*traxels_tracklet_map)[tn_next].begin(),
                                           (*traxels_tracklet_map)[tn_next].end());
                        }

                        for (std::vector<Traxel>::const_iterator tr = traxels.begin(); tr != traxels.end(); ++tr) {
                            timesteps.push_back(tr->Timestep);
                            tracklet.push_back(*tr);
                        }
                        active_outgoing = true;
                    }
                }
                if (active_outgoing) {
                    tn = tn_next;
                } else { // no active outgoing arc found -- tracklet ends
                    break;
                }
            }

            HypothesesGraph::Node curr_tracklet_node = tracklet_graph.add_node(timesteps);
            tracklet_map.set(curr_tracklet_node, tracklet);

            // store the outgoing arcs of the traxel node
            if (tn_outarcs.size() > 0) {
                for(std::vector<HypothesesGraph::Arc>::const_iterator arc_it = tn_outarcs.begin(); arc_it!=tn_outarcs.end(); ++arc_it) {
                    HypothesesGraph::Node traxel_node_to = traxel_graph.target(*arc_it);
                    if (node_link_map.find(traxel_node_to) == node_link_map.end()) {
                        node_link_map.insert(
                                    std::pair<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > (
                                        traxel_node_to, std::vector<HypothesesGraph::Node>() ) );
                    }
                    node_link_map[traxel_node_to].push_back(curr_tracklet_node);
                }
            }

            // set the incoming arcs of the tracklet (traxel_node is the first node in the tracklet)
            if (!traxel_nodes_are_tracklets) {
                assert(traxel_map[traxel_node].Id==tracklet[0].Id);
                assert(traxel_map[traxel_node].Timestep==tracklet[0].Timestep);
            } else {
                assert((*traxels_tracklet_map)[traxel_node][0].Id==tracklet[0].Id);
                assert((*traxels_tracklet_map)[traxel_node][0].Timestep==tracklet[0].Timestep);
            }
            std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> >::const_iterator traxel_to_it = node_link_map.find(traxel_node);

            if (traxel_to_it == node_link_map.end()) {
                LOG(logDEBUG3) << "There is no incoming arc for traxel_node " << traxel_graph.id(traxel_node);
            } else {
                for(std::vector<HypothesesGraph::Node>::const_iterator tracklet_from_it = (*traxel_to_it).second.begin();
                    tracklet_from_it != (*traxel_to_it).second.end(); ++tracklet_from_it) {
                    tracklet_graph.addArc(*tracklet_from_it,curr_tracklet_node);
                }
            }
        }
    }
}

namespace {
std::vector<HypothesesGraph::Arc> getOutgoingArcs(const HypothesesGraph& graph, const HypothesesGraph::Node& n) {
    std::vector<HypothesesGraph::Arc> result;
    for(HypothesesGraph::OutArcIt a(graph, n); a != lemon::INVALID; ++a) {
        result.push_back(a);
    }
    return result;
}

std::vector<HypothesesGraph::Arc> getIncomingArcs(const HypothesesGraph& graph, const HypothesesGraph::Node& n) {
    std::vector<HypothesesGraph::Arc> result;
    for(HypothesesGraph::InArcIt a(graph, n); a != lemon::INVALID; ++a) {
        result.push_back(a);
    }
    return result;
}

void addNodeToTracklet(const HypothesesGraph& traxel_graph, HypothesesGraph& tracklet_graph,
                       const HypothesesGraph::Node& traxel_node, const HypothesesGraph::Node& ancestor_traxel_node,
                       std::map<HypothesesGraph::Node, HypothesesGraph::Node>& traxel2tracklet, double traxel_arc_dist,
                       std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> >& tracklet2traxel, const int arc_id) {
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = traxel_graph.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = tracklet_graph.get(node_tracklet());
    property_map<tracklet_intern_dist, HypothesesGraph::base_graph>::type& tracklet_arc_dist_map = tracklet_graph.get(tracklet_intern_dist());
    property_map<tracklet_intern_arc_ids, HypothesesGraph::base_graph>::type& tracklet_arc_id_map = tracklet_graph.get(tracklet_intern_arc_ids());

    //	property_map<node_timestep, HypothesesGraph::base_graph>::type& timestep_map = tracklet_graph.get(node_timestep());

    assert(traxel2tracklet.find(ancestor_traxel_node) != traxel2tracklet.end());
    HypothesesGraph::Node tracklet_node = traxel2tracklet[ancestor_traxel_node];
    std::vector<Traxel> tracklet = tracklet_map[tracklet_node];
    Traxel tr = traxel_map[traxel_node];
    tracklet.push_back(tr);
    tracklet_map.set(tracklet_node, tracklet);
    traxel2tracklet[traxel_node] = tracklet_node;


  	if(traxel_graph.getProperties().count("appearance_label") > 0){
		// create appropriate labels in tracklet node
			tracklet_graph.add_appearance_label(tracklet_node,traxel_graph.get(appearance_label())[traxel_node]);
			tracklet_graph.add_disappearance_label(tracklet_node,traxel_graph.get(disappearance_label())[traxel_node]);
			tracklet_graph.add_division_label(tracklet_node,traxel_graph.get(division_label())[traxel_node]);
	}
    //	size_t timestep = tr.Timestep;
    //	timestep_map[tracklet_node].add(timestep);
    //	tracklet_graph.timesteps_.insert(timestep);

    // add internal arc
    std::vector<double> arc_dists = tracklet_arc_dist_map[tracklet_node];
    arc_dists.push_back(traxel_arc_dist);
    tracklet_arc_dist_map.set(tracklet_node,arc_dists);

    tracklet2traxel[tracklet_node].push_back(traxel_node);

    std::vector<int> arc_ids = tracklet_arc_id_map[tracklet_node];
    arc_ids.push_back(arc_id);
    tracklet_arc_id_map.set(tracklet_node, arc_ids);
    LOG(logDEBUG4) << "addNodeToTracklet: added arc_id " << arc_id;
}

void addNodeToGraph(const HypothesesGraph& traxel_graph, HypothesesGraph& tracklet_graph,
                    const HypothesesGraph::Node& traxel_node, std::map<HypothesesGraph::Node, HypothesesGraph::Node>& traxel2tracklet,
                    std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> >& tracklet2traxel,
                    bool with_ground_truth_labeling=false) {
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = traxel_graph.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = tracklet_graph.get(node_tracklet());
    property_map<tracklet_intern_dist, HypothesesGraph::base_graph>::type& tracklet_intern_dist_map = tracklet_graph.get(tracklet_intern_dist());
    std::vector<Traxel> tracklet;
    Traxel tr = traxel_map[traxel_node];
    tracklet.push_back(tr);
    size_t timestep = tr.Timestep;
    HypothesesGraph::Node tracklet_node = tracklet_graph.add_node(timestep);
    LOG(logDEBUG4) << "added tracklet node " << tracklet_graph.id(tracklet_node);
    tracklet_map.set(tracklet_node, tracklet);
    traxel2tracklet[traxel_node] = tracklet_node;
    std::vector<double> arc_dists;
    tracklet_intern_dist_map.set(tracklet_node, arc_dists);

    if(with_ground_truth_labeling){
		// create appropriate labels in tracklet node
		tracklet_graph.add_appearance_label(tracklet_node,traxel_graph.get(appearance_label())[traxel_node]);
		tracklet_graph.add_disappearance_label(tracklet_node,traxel_graph.get(disappearance_label())[traxel_node]);
		tracklet_graph.add_division_label(tracklet_node,traxel_graph.get(division_label())[traxel_node]);
	}

    tracklet2traxel[tracklet_node].push_back(traxel_node);
}

void addArcsToGraph(const HypothesesGraph& traxel_graph, HypothesesGraph& tracklet_graph,
                    const std::vector<HypothesesGraph::Arc>& incoming_arcs, std::map<HypothesesGraph::Node, HypothesesGraph::Node>& traxel2tracklet,
                    bool with_ground_truth_labeling=false) {
    property_map<arc_distance, HypothesesGraph::base_graph>::type& traxel_arc_distances = traxel_graph.get(arc_distance());
    property_map<arc_distance, HypothesesGraph::base_graph>::type& tracklet_arc_distances = tracklet_graph.get(arc_distance());
    property_map<traxel_arc_id, HypothesesGraph::base_graph>::type& traxel_arc_ids = tracklet_graph.get(traxel_arc_id());

    for(std::vector<HypothesesGraph::Arc>::const_iterator a = incoming_arcs.begin(); a!=incoming_arcs.end(); ++a) {
        HypothesesGraph::Arc arc = *a;
        LOG(logDEBUG4) << "traxel arc source " << traxel_graph.id(traxel_graph.source(arc));
        LOG(logDEBUG4) << "traxel arc target " << traxel_graph.id(traxel_graph.target(arc));
        HypothesesGraph::Node from = traxel2tracklet[traxel_graph.source(arc)];
        HypothesesGraph::Node to = traxel2tracklet[traxel_graph.target(arc)];
        LOG(logDEBUG4) << "tracklet node from " << tracklet_graph.id(traxel2tracklet[traxel_graph.source(arc)]);
        LOG(logDEBUG4) << "tracklet node to " << tracklet_graph.id(traxel2tracklet[traxel_graph.target(arc)]);
        assert(from != to);

        HypothesesGraph::Arc tracklet_arc = tracklet_graph.addArc(from, to);
        if(with_ground_truth_labeling){
			tracklet_graph.add_arc_label(tracklet_arc,traxel_graph.get(arc_label())[arc]);
		}
        tracklet_arc_distances.set(tracklet_arc,traxel_arc_distances[arc]);
        traxel_arc_ids.set(tracklet_arc, (int) traxel_graph.id(arc));
    }
}

} // anonymous namespace

std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > generateTrackletGraph2(const HypothesesGraph& traxel_graph, HypothesesGraph& tracklet_graph) {
    property_map<arc_distance, HypothesesGraph::base_graph>::type& traxel_arc_dist_map = traxel_graph.get(arc_distance());

    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    node_timestep_map_t& node_timestep_map = traxel_graph.get(node_timestep());

    // add empty traxel_map to the tracklet graph in order to make the tracklet graph equivalent to traxelgraphs
    tracklet_graph.add(node_traxel()).add(arc_distance());

    tracklet_graph.add(node_tracklet()).add(tracklet_intern_dist()).add(tracklet_intern_arc_ids()).add(traxel_arc_id());

    std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > tracklet_node_to_traxel_nodes;
    // maps traxel_nodes to tracklet_nodes
    std::map<HypothesesGraph::Node, HypothesesGraph::Node > traxel_node_to_tracklet_node;

    bool with_ground_truth_labeling = (traxel_graph.getProperties().count("appearance_label") > 0);
    if (with_ground_truth_labeling) {
        tracklet_graph.add(appearance_label()).add(disappearance_label()).add(division_label()).add(arc_label());
    }

    for(int t = traxel_graph.earliest_timestep(); t <= traxel_graph.latest_timestep(); ++t) {
        for(node_timestep_map_t::ItemIt traxel_node(node_timestep_map, t); traxel_node!=lemon::INVALID; ++traxel_node) {
            LOG(logDEBUG4) << "traxel_node = " << traxel_graph.id(traxel_node);
            vector<HypothesesGraph::Arc> incoming_arcs = getIncomingArcs(traxel_graph, traxel_node);
            if(incoming_arcs.size() != 1) {
                addNodeToGraph(traxel_graph, tracklet_graph, traxel_node, traxel_node_to_tracklet_node, tracklet_node_to_traxel_nodes,
                               with_ground_truth_labeling);
                addArcsToGraph(traxel_graph, tracklet_graph, incoming_arcs, traxel_node_to_tracklet_node, with_ground_truth_labeling);
                LOG(logDEBUG4) << "traxel2tracklet.size(): " << traxel_node_to_tracklet_node.size();
                continue;
            }

            HypothesesGraph::Node ancestor = traxel_graph.source(incoming_arcs[0]);
            if(getOutgoingArcs(traxel_graph, ancestor).size() > 1) {
                addNodeToGraph(traxel_graph, tracklet_graph, traxel_node, traxel_node_to_tracklet_node, tracklet_node_to_traxel_nodes,
                               with_ground_truth_labeling);
                assert(incoming_arcs.size() == 1);
                addArcsToGraph(traxel_graph, tracklet_graph, incoming_arcs, traxel_node_to_tracklet_node,
                               with_ground_truth_labeling);
                LOG(logDEBUG4) << "traxel2tracklet.size(): " << traxel_node_to_tracklet_node.size();
                continue;
            }

            double dist = traxel_arc_dist_map[incoming_arcs[0]];
            int arc_id = traxel_graph.id(incoming_arcs[0]);
            addNodeToTracklet(traxel_graph, tracklet_graph, traxel_node, ancestor, traxel_node_to_tracklet_node, dist,
                              tracklet_node_to_traxel_nodes, arc_id);
        }
    }

    return tracklet_node_to_traxel_nodes;
}




//
// state_of_nodes()
//
boost::shared_ptr<std::vector< std::map<unsigned int, bool> > > state_of_nodes(const HypothesesGraph& g) {
    LOG(logDEBUG) << "detections(): entered";
    boost::shared_ptr<vector< map<unsigned int, bool> > > ret(new vector< map<unsigned int, bool> >);

    // required node properties: timestep, traxel, active
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    node_timestep_map_t& node_timestep_map = g.get(node_timestep());
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map_t;
    node_traxel_map_t& node_traxel_map = g.get(node_traxel());
    property_map<node_active, HypothesesGraph::base_graph>::type* node_active_map;
    property_map<node_active2, HypothesesGraph::base_graph>::type* node_active2_map;
    bool active2_used = false;
    if (g.getProperties().count("node_active") > 0) {
        node_active_map = &g.get(node_active());
    } else if (g.getProperties().count("node_active2") > 0) {
        node_active2_map = &g.get(node_active2());
        active2_used = true;
    }

    // for every timestep
    for(int t = g.earliest_timestep(); t <= g.latest_timestep(); ++t) {
        ret->push_back(map<unsigned int, bool>());
        for(node_timestep_map_t::ItemIt node_at(node_timestep_map, t); node_at!=lemon::INVALID; ++node_at) {
            assert(node_traxel_map[node_at].Timestep == t);
            unsigned int id = node_traxel_map[node_at].Id;
            bool active = false;
            if (active2_used) {
                active = ((*node_active2_map)[node_at] > 0);
            } else {
                active = (*node_active_map)[node_at];
            }
            (*ret)[t-g.earliest_timestep()][id] = active;
        }
    }

    return ret;
}


//
// write_lgf()
//

namespace {
template <typename T>
struct TypeToStrConverter {
	std::string operator()(const T& t) {
		stringstream ss;
		boost::archive::text_oarchive oa(ss);
		oa & t;
		return ss.str();
	}
};
}
  void write_lgf(const HypothesesGraph& g, std::ostream& os, std::map<std::string, bool>& config ) {
    lemon::DigraphWriter<HypothesesGraph> writer( g, os );    

    if(config["node_timestep"])
        writer.nodeMap("node_timestep", g.get(node_timestep()));

    if(config["node_active"])
        writer.nodeMap("node_active", g.get(node_active()));

    if(config["node_active2"])
        writer.nodeMap("node_active2", g.get(node_active2()));

    if(config["node_active_count"])
        writer.nodeMap("node_active_count", g.get(node_active_count()), TypeToStrConverter<vector<long unsigned int> >());

    if(config["node_offered"])
        writer.nodeMap("node_offered", g.get(node_offered()));

    if(config["split_from"])
        writer.nodeMap("split_from", g.get(split_from()));

    if(config["division_active"])
        writer.nodeMap("division_active", g.get(division_active()));

    if(config["merger_resolved_to"])
        writer.nodeMap("merger_resolved_to", g.get(merger_resolved_to()), TypeToStrConverter<vector<unsigned int> >());

    if(config["node_originated_from"])
        writer.nodeMap("node_originated_from", g.get(node_originated_from()), TypeToStrConverter<vector<unsigned int> >());

    if(config["node_resolution_candidate"])
        writer.nodeMap("node_resolution_candidate", g.get(node_resolution_candidate()));

    if(config["arc_distance"])
        writer.arcMap("arc_distance", g.get(arc_distance()));

    if(config["traxel_arc_id"])
        writer.arcMap("traxel_arc_id", g.get(traxel_arc_id()));

    if(config["arc_vol_ratio"])
        writer.arcMap("arc_vol_ratio", g.get(arc_vol_ratio()));

    if(config["arc_from_timestep"])
        writer.arcMap("arc_from_timestep", g.get(arc_from_timestep()));

    if(config["arc_to_timestep"])
        writer.arcMap("arc_to_timestep", g.get(arc_to_timestep()));

    if(config["arc_active"])
        writer.arcMap("arc_active", g.get(arc_active()));

    if(config["arc_resolution_candidate"])
        writer.arcMap("arc_resolution_candidate", g.get(arc_resolution_candidate()));

    if(config["tracklet_intern_dist"])
        writer.nodeMap("tracklet_intern_dist", g.get(tracklet_intern_dist()), TypeToStrConverter<vector<double> >());

    if(config["tracklet_intern_arc_ids"])
        writer.nodeMap("tracklet_intern_arc_ids", g.get(tracklet_intern_arc_ids()), TypeToStrConverter<vector<int> >());

    if(config["arc_active_count"])
        writer.arcMap("arc_active_count", g.get(arc_active_count()), TypeToStrConverter<vector<bool> >());

    if(config["node_traxel"])
    	writer.nodeMap("node_traxel", g.get(node_traxel()), TypeToStrConverter<Traxel>());

    if(config["node_tracklet"])
        writer.nodeMap("node_tracklet", g.get(node_tracklet()), TypeToStrConverter<vector<Traxel>>());

    writer.run();
  }
  //
  // read_lgf()
  //
  namespace {
  template <typename T>
  struct StrToTypeConverter {
	  T operator()(const string& s) {
		  stringstream ss(s);
		  boost::archive::text_iarchive ia(ss);
		  T t;
		  ia & t;
		  return t;
	  }
  };
  }

  void read_lgf( HypothesesGraph& g, std::istream& is, std::map<std::string, bool>& config) {
	  lemon::DigraphReader<HypothesesGraph> reader( g, is );

    if(config["node_timestep"])
    {
        g.add(node_timestep());
        reader.nodeMap("node_timestep", g.get(node_timestep()));
    }

    if(config["node_active"])
    {
        g.add(node_active());
        reader.nodeMap("node_active", g.get(node_active()));
    }

    if(config["node_active2"])
    {
        g.add(node_active2());
        reader.nodeMap("node_active2", g.get(node_active2()));
    }

    if(config["node_active_count"])
    {
        g.add(node_active_count());
        reader.nodeMap("node_active_count", g.get(node_active_count()), StrToTypeConverter<vector<long unsigned int> >());
    }

    if(config["node_offered"])
    {
        g.add(node_offered());
        reader.nodeMap("node_offered", g.get(node_offered()));
    }

    if(config["split_from"])
    {
        g.add(split_from());
        reader.nodeMap("split_from", g.get(split_from()));
    }

    if(config["division_active"])
    {
        g.add(division_active());
        reader.nodeMap("division_active", g.get(division_active()));
    }

    if(config["merger_resolved_to"])
    {
        g.add(merger_resolved_to());
        reader.nodeMap("merger_resolved_to", g.get(merger_resolved_to()), StrToTypeConverter<vector<unsigned int> >());
    }

    if(config["node_originated_from"])
    {
        g.add(node_originated_from());
        reader.nodeMap("node_originated_from", g.get(node_originated_from()), StrToTypeConverter<vector<unsigned int> >());
    }

    if(config["node_resolution_candidate"])
    {
        g.add(node_resolution_candidate());
        reader.nodeMap("node_resolution_candidate", g.get(node_resolution_candidate()));
    }

    if(config["arc_distance"])
    {
        g.add(arc_distance());
        reader.arcMap("arc_distance", g.get(arc_distance()));
    }

    if(config["traxel_arc_id"])
    {
        g.add(traxel_arc_id());
        reader.arcMap("traxel_arc_id", g.get(traxel_arc_id()));
    }

    if(config["arc_vol_ratio"])
    {
        g.add(arc_vol_ratio());
        reader.arcMap("arc_vol_ratio", g.get(arc_vol_ratio()));
    }

    if(config["arc_from_timestep"])
    {
        g.add(arc_from_timestep());
        reader.arcMap("arc_from_timestep", g.get(arc_from_timestep()));
    }

    if(config["arc_to_timestep"])
    {
        g.add(arc_to_timestep());
        reader.arcMap("arc_to_timestep", g.get(arc_to_timestep()));
    }

    if(config["arc_active"])
    {
        g.add(arc_active());
        reader.arcMap("arc_active", g.get(arc_active()));
    }

    if(config["arc_resolution_candidate"])
    {
        g.add(arc_resolution_candidate());
        reader.arcMap("arc_resolution_candidate", g.get(arc_resolution_candidate()));
    }

    if(config["tracklet_intern_dist"])
    {
        g.add(tracklet_intern_dist());
        reader.nodeMap("tracklet_intern_dist", g.get(tracklet_intern_dist()), StrToTypeConverter<vector<double> >());
    }

    if(config["tracklet_intern_arc_ids"])
    {
        g.add(tracklet_intern_arc_ids());
        reader.nodeMap("tracklet_intern_arc_ids", g.get(tracklet_intern_arc_ids()), StrToTypeConverter<vector<int> >());
    }

    if(config["arc_active_count"])
    {
    	g.add(arc_active_count());
    	reader.arcMap("arc_active_count", g.get(arc_active_count()), StrToTypeConverter<vector<bool> >());
    }

    if(config["node_traxel"])
    {
    	g.add(node_traxel());
    	reader.nodeMap("node_traxel", g.get(node_traxel()), StrToTypeConverter<Traxel>());
    }

    if(config["node_tracklet"])
    {
        g.add(node_tracklet());
        reader.nodeMap("node_tracklet", g.get(node_tracklet()), StrToTypeConverter<vector<Traxel>>());
    }

    reader.run();
  }



////
//// class HypothesesBuilder
////
HypothesesGraph* HypothesesBuilder::build() const {
    // construct an empty HypothesesGraph with all desired additional
    // properties added
    HypothesesGraph* graph = construct();
    // add object nodes and set node properties
    graph = add_nodes(graph);
    // connect object nodes and set edge properties
    graph = add_edges(graph);

    return graph;
}


namespace {
double getDivisionProbability(const Traxel& tr) {
    FeatureMap::const_iterator it = tr.features.find("divProb");
    if (it == tr.features.end()) {
        throw runtime_error("getDivisionProbability(): divProb feature not in traxel");
    }
    return it->second[0];
}
}

////
//// class SingleTimestepTraxel_HypothesesBuilder
////
HypothesesGraph* SingleTimestepTraxel_HypothesesBuilder::construct() const {
    HypothesesGraph* graph = new HypothesesGraph();
    // store traxels inside the graph data structure
    graph->add(node_traxel());
    return graph;
}

HypothesesGraph* SingleTimestepTraxel_HypothesesBuilder::add_nodes(HypothesesGraph* graph) const {
    LOG(logDEBUG) << "SingleTimestepTraxel_HypothesesBuilder::add_nodes(): entered";
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_m = graph->get(node_traxel());

    for(TraxelStoreByTimestep::const_iterator it = ts_->begin(); it!= ts_->end(); ++it) {
        HypothesesGraph::Node node = graph->add_node(it->Timestep);
        traxel_m.set(node, *it);
    }

    return graph;
}

HypothesesGraph* SingleTimestepTraxel_HypothesesBuilder::add_edges(
        HypothesesGraph* graph) const {
    LOG(logDEBUG) << "SingleTimestepTraxel_HypothesesBuilder::add_edges(): entered";
    typedef HypothesesGraph::node_timestep_map::Value timestep_t;
    const set<timestep_t>& timesteps = graph->timesteps();
    // iterate over all timesteps except the last
    for (set<timestep_t>::const_iterator t = timesteps.begin();
         t != (--timesteps.end()); ++t) {
        add_edges_at(graph, *t);
    }

    // if the forward_backward option is enabled, go again through the graph, this time
    // reversely and add the nearest neighbors if not already present
    if (options_.forward_backward) {
        // reversely iterate over all timesteps except the first
        for (set<timestep_t>::const_reverse_iterator t = timesteps.rbegin();
             t != (--timesteps.rend()); ++t) {
            add_edges_at(graph, *t, true);
        }
    }

    return graph;
}

HypothesesGraph* SingleTimestepTraxel_HypothesesBuilder::add_edges_at(HypothesesGraph* graph,
                                                                      int timestep, bool reverse) const {
    const HypothesesGraph::node_timestep_map& timemap = graph->get(
                node_timestep());
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type traxelmap_t;
    const traxelmap_t& traxelmap = graph->get(node_traxel());
    const TraxelStoreByTimeid& traxels_by_timeid = ts_->get<by_timeid>();
    const TraxelStoreByTimestep& traxels_by_timestep = ts_->get<by_timestep>();

    int to_timestep = timestep + 1;
    if (reverse) {
        // iterating through the graph backward in time
        to_timestep = timestep - 1;
    }

    //// find k nearest neighbors in next timestep
    // init nearest neighbor search
    pair<TraxelStoreByTimestep::const_iterator,
            TraxelStoreByTimestep::const_iterator> traxels_at =
            traxels_by_timestep.equal_range(to_timestep);

    for (TraxelStoreByTimestep::const_iterator it = traxels_at.first;
         it != traxels_at.second; ++it) {
        assert(it->Timestep == to_timestep);
    }

    NearestNeighborSearch nns(traxels_at.first, traxels_at.second, reverse);


    // establish transition edges between a current node and appropriate nodes in next timestep
    for (HypothesesGraph::node_timestep_map::ItemIt curr_node(timemap,
                                                              timestep); curr_node != lemon::INVALID; ++curr_node) {
        assert(timemap[curr_node] == timestep);
        assert(traxelmap[curr_node].Timestep == timestep);

        // if we want to consider divisions already in the Hypotheses graph
        // make sure that each potentially dividing cell has 2 nearest neighbors
        // (but only if we go through the graph forward in time)
        unsigned int max_nn = options_.max_nearest_neighbors;
        if (options_.consider_divisions && !reverse && max_nn < 2) {
            double div_prob = getDivisionProbability(traxelmap[curr_node]);
            if (div_prob > options_.division_threshold) {
                max_nn = 2;
            }
        }

        // search
        map<unsigned int, double> nearest_neighbors = nns.knn_in_range(
                    traxelmap[curr_node], options_.distance_threshold,
                    max_nn, reverse);

        //// connect current node with k nearest neighbor nodes
        for (map<unsigned int, double>::const_iterator neighbor =
             nearest_neighbors.begin(); neighbor != nearest_neighbors.end();
             ++neighbor) {
            // connect with one of the neighbor nodes
            TraxelStoreByTimeid::iterator neighbor_traxel =
                    traxels_by_timeid.find(
                        boost::make_tuple(to_timestep, neighbor->first));
            assert(neighbor_traxel->Timestep == to_timestep);
            assert(neighbor_traxel->Timestep != traxelmap[curr_node].Timestep);
            traxelmap_t::ItemIt neighbor_node(traxelmap, *neighbor_traxel);
            assert(curr_node != neighbor_node);
            if (!reverse) {
                // if we go through the graph forward in time, add an arc from curr_node to neighbor_node
                LOG(logDEBUG4) << "added arc from traxel " << traxelmap[curr_node].Id << " to " <<
                                  traxelmap[neighbor_node].Id;
                graph->addArc(curr_node, neighbor_node);
            } else {
                // if we go through the graph backward in time, add an arc from neighbor_node to curr_node
                // if not already present
                if (lemon::findArc(*graph,neighbor_node,curr_node) == lemon::INVALID) {
                    graph->addArc(neighbor_node, curr_node);
                    LOG(logDEBUG4) << "added backward arc from traxel " << traxelmap[neighbor_node].Id << " to " <<
                                      traxelmap[curr_node].Id;
                }
            }
        }
    }

    return graph;
}

// graph copy methods
template<class Graph>
void PropertyGraph<Graph>::copy(PropertyGraph<Graph>& src, PropertyGraph<Graph>& dest)
{
    throw std::runtime_error("Can only copy lemon::ListDigraph");
}

template<>
void PropertyGraph<lemon::ListDigraph>::copy(PropertyGraph<lemon::ListDigraph>& src, PropertyGraph<lemon::ListDigraph>& dest)
{
    lemon::DigraphCopy<PropertyGraph<lemon::ListDigraph>, PropertyGraph<lemon::ListDigraph> > graph_copy(src, dest);

    // create copies of all property maps
    if(src.has_property(node_timestep()))
    {
        dest.add(node_timestep());
        graph_copy.nodeMap(src.get(node_timestep()), dest.get(node_timestep()));
    }

    if(src.has_property(node_active()))
    {
        dest.add(node_active());
        graph_copy.nodeMap(src.get(node_active()), dest.get(node_active()));
    }

    if(src.has_property(node_active2()))
    {
        dest.add(node_active2());
        graph_copy.nodeMap(src.get(node_active2()), dest.get(node_active2()));
    }

    if(src.has_property(node_active_count()))
    {
        dest.add(node_active_count());
        graph_copy.nodeMap(src.get(node_active_count()), dest.get(node_active_count()));
    }

    if(src.has_property(node_offered()))
    {
        dest.add(node_offered());
        graph_copy.nodeMap(src.get(node_offered()), dest.get(node_offered()));
    }

    if(src.has_property(split_from()))
    {
        dest.add(split_from());
        graph_copy.nodeMap(src.get(split_from()), dest.get(split_from()));
    }

    if(src.has_property(division_active()))
    {
        dest.add(division_active());
        graph_copy.nodeMap(src.get(division_active()), dest.get(division_active()));
    }

    if(src.has_property(merger_resolved_to()))
    {
        dest.add(merger_resolved_to());
        graph_copy.nodeMap(src.get(merger_resolved_to()), dest.get(merger_resolved_to()));
    }

    if(src.has_property(node_originated_from()))
    {
        dest.add(node_originated_from());
        graph_copy.nodeMap(src.get(node_originated_from()), dest.get(node_originated_from()));
    }

    if(src.has_property(node_resolution_candidate()))
    {
        dest.add(node_resolution_candidate());
        graph_copy.nodeMap(src.get(node_resolution_candidate()), dest.get(node_resolution_candidate()));
    }

    if(src.has_property(arc_distance()))
    {
        dest.add(arc_distance());
        graph_copy.arcMap(src.get(arc_distance()), dest.get(arc_distance()));
    }

    if(src.has_property(traxel_arc_id()))
    {
        dest.add(traxel_arc_id());
        graph_copy.arcMap(src.get(traxel_arc_id()), dest.get(traxel_arc_id()));
    }

    if(src.has_property(arc_vol_ratio()))
    {
        dest.add(arc_vol_ratio());
        graph_copy.arcMap(src.get(arc_vol_ratio()), dest.get(arc_vol_ratio()));
    }

    if(src.has_property(arc_from_timestep()))
    {
        dest.add(arc_from_timestep());
        graph_copy.arcMap(src.get(arc_from_timestep()), dest.get(arc_from_timestep()));
    }

    if(src.has_property(arc_to_timestep()))
    {
        dest.add(arc_to_timestep());
        graph_copy.arcMap(src.get(arc_to_timestep()), dest.get(arc_to_timestep()));
    }

    if(src.has_property(arc_active()))
    {
        dest.add(arc_active());
        graph_copy.arcMap(src.get(arc_active()), dest.get(arc_active()));
    }

    if(src.has_property(arc_resolution_candidate()))
    {
        dest.add(arc_resolution_candidate());
        graph_copy.arcMap(src.get(arc_resolution_candidate()), dest.get(arc_resolution_candidate()));
    }

    if(src.has_property(tracklet_intern_dist()))
    {
        dest.add(tracklet_intern_dist());
        graph_copy.nodeMap(src.get(tracklet_intern_dist()), dest.get(tracklet_intern_dist()));
    }

    if(src.has_property(tracklet_intern_arc_ids()))
    {
        dest.add(tracklet_intern_arc_ids());
        graph_copy.nodeMap(src.get(tracklet_intern_arc_ids()), dest.get(tracklet_intern_arc_ids()));
    }

    if(src.has_property(arc_active_count()))
    {
        dest.add(arc_active_count());
        graph_copy.arcMap(src.get(arc_active_count()), dest.get(arc_active_count()));
    }

    if(src.has_property(node_traxel()))
    {
        dest.add(node_traxel());
        graph_copy.nodeMap(src.get(node_traxel()), dest.get(node_traxel()));
    }

    graph_copy.run();
}

void HypothesesGraph::copy(HypothesesGraph &src, HypothesesGraph &dest)
{
    // copy "parents"
    PropertyGraph<lemon::ListDigraph>::copy(src, dest);

    dest.timesteps_ = src.timesteps_;
}

} /* namespace pgmlink */

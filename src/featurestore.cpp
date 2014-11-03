#include "pgmlink/featurestore.h"
#include "pgmlink/traxels.h"
#include <iostream>

namespace pgmlink {

FeatureStore::FeatureStore()
{
    std::cout << "FeatureStore created" << std::endl;
}

FeatureMap &FeatureStore::get_traxel_features(int timestep, unsigned int id)
{
    return traxel_feature_map_[{std::make_pair(timestep, id)}];
}

FeatureMap &FeatureStore::get_traxel_features(const Traxel &traxel)
{
    return traxel_feature_map_[{std::make_pair(traxel.Timestep, traxel.Id)}];
}

FeatureMap &FeatureStore::get_traxel_features(const Traxel &traxel_a, const Traxel &traxel_b)
{
    return traxel_feature_map_[{std::make_pair(traxel_a.Timestep, traxel_a.Id),
            std::make_pair(traxel_b.Timestep, traxel_b.Id)}];
}

FeatureMap &FeatureStore::get_traxel_features(const Traxel &traxel_a, const Traxel &traxel_b, const Traxel &traxel_c)
{
    return traxel_feature_map_[{std::make_pair(traxel_a.Timestep, traxel_a.Id),
            std::make_pair(traxel_b.Timestep, traxel_b.Id),
            std::make_pair(traxel_c.Timestep, traxel_c.Id)}];
}

FeatureMap &FeatureStore::get_traxel_features(const std::vector<const Traxel *> &traxels)
{
    std::vector<TimeId> keys;

    for(const auto t : traxels)
    {
        keys.push_back(std::make_pair(t->Timestep, t->Id));
    }
}

void FeatureStore::dump(std::ostream& stream)
{
    for(TraxelFeatureMap::iterator it = traxel_feature_map_.begin(); it != traxel_feature_map_.end(); ++it)
    {
        stream << "Traxel (" << it->first[0].first << ", " << it->first[0].second << ")\n";
        FeatureMap& feature_map = it->second;
        for(FeatureMap::iterator map_it = feature_map.begin(); map_it != feature_map.end(); ++map_it)
        {
            stream << "\tFeature \"" << map_it->first << "\": ";
            feature_array& feat = map_it->second;
            for(auto f : feat)
            {
                stream << f << " ";
            }
            stream << "\n";
        }
    }
    stream << std::endl;
}

void FeatureStore::dump(int timestep, unsigned int id, std::ostream &stream)
{
    TraxelFeatureMap::iterator it = traxel_feature_map_.find({std::make_pair(timestep, id)});

    stream << "Traxel (" << it->first[0].first << ", " << it->first[0].second << ")\n";
    FeatureMap& feature_map = it->second;
    for(FeatureMap::iterator map_it = feature_map.begin(); map_it != feature_map.end(); ++map_it)
    {
        stream << "\tFeature \"" << map_it->first << "\": ";
        feature_array& feat = map_it->second;
        for(auto f : feat)
        {
            stream << f << " ";
        }
        stream << "\n";
    }
}

} // end namespace pgmlink

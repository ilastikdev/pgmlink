/**
   @file
   @ingroup tracking
   @brief tracking datastructures
*/

#ifndef TRAXELS_H
#define TRAXELS_H

#include <set>
#include <string>
#include <iostream>
#include <ostream>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "pgmlink/pgmlink_export.h"
#include "pgmlink/log.h"
#include "featurestore.h"

namespace pgmlink {

/*
 * Feature datastructures now in featurestore.h
typedef float feature_type;
typedef std::vector<feature_type> feature_array;
typedef std::vector<feature_array> feature_arrays;
typedef std::map<std::string,feature_array> FeatureMap;
*/

//
// retrieve spatial coordinates from features
//
class Locator {
 public:
  PGMLINK_EXPORT Locator( std::string fn,
	                      double x_scale = 1.0,
	                      double y_scale = 1.0,
	                      double z_scale = 1.0 ) 
  : x_scale(x_scale), y_scale(y_scale), z_scale(z_scale), feature_name_(fn) 
  {}
 
  PGMLINK_EXPORT virtual Locator* clone() = 0;
  PGMLINK_EXPORT virtual ~Locator() {};

  PGMLINK_EXPORT virtual bool is_applicable(const FeatureMap& m) const { return m.count(feature_name_)==1; };
  PGMLINK_EXPORT virtual double X(const FeatureMap&) const = 0;
  PGMLINK_EXPORT virtual double Y(const FeatureMap&) const = 0;
  PGMLINK_EXPORT virtual double Z(const FeatureMap&) const = 0;

  double x_scale, y_scale, z_scale;

 protected:
  std::string feature_name_;
  PGMLINK_EXPORT double coordinate_from(const FeatureMap&, size_t idx) const;

 private:
  // boost serialize
  friend class boost::serialization::access;
  template< typename Archive >
    void serialize( Archive&, const unsigned int /*version*/ );
};

class ComLocator 
: public Locator 
{
 public:
  PGMLINK_EXPORT ComLocator() 
  : Locator("com")
  {}
  
  PGMLINK_EXPORT virtual ComLocator* clone() { return new ComLocator(*this); }
  PGMLINK_EXPORT double X(const FeatureMap& m) const { return x_scale * coordinate_from(m, 0); }
  PGMLINK_EXPORT double Y(const FeatureMap& m) const { return y_scale * coordinate_from(m, 1); }
  PGMLINK_EXPORT double Z(const FeatureMap& m) const { return z_scale * coordinate_from(m, 2); }

 private:
  // boost serialize
  friend class boost::serialization::access;
  template< typename Archive >
    void serialize( Archive&, const unsigned int /*version*/ );
};

class ComCorrLocator 
: public Locator 
{
 public:
  PGMLINK_EXPORT ComCorrLocator() 
  : Locator("com_corrected") 
  {}

  PGMLINK_EXPORT virtual ComCorrLocator* clone() { return new ComCorrLocator(*this); }
  PGMLINK_EXPORT double X(const FeatureMap& m) const { return x_scale * coordinate_from(m, 0); }
  PGMLINK_EXPORT double Y(const FeatureMap& m) const { return y_scale * coordinate_from(m, 1); }
  PGMLINK_EXPORT double Z(const FeatureMap& m) const { return z_scale * coordinate_from(m, 2); }

 private:
  // boost serialize
  friend class boost::serialization::access;
  template< typename Archive >
    void serialize( Archive&, const unsigned int /*version*/ );
};

class IntmaxposLocator 
: public Locator 
{
 public:
  PGMLINK_EXPORT IntmaxposLocator() 
  : Locator("intmaxpos") 
  {}

  PGMLINK_EXPORT virtual IntmaxposLocator* clone() { return new IntmaxposLocator(*this); }
  PGMLINK_EXPORT double X(const FeatureMap& m) const { return x_scale * coordinate_from(m, 1); }
  PGMLINK_EXPORT double Y(const FeatureMap& m) const { return y_scale * coordinate_from(m, 2); }
  PGMLINK_EXPORT double Z(const FeatureMap& m) const { return z_scale * coordinate_from(m, 3); }

 private:
  // boost serialize
  friend class boost::serialization::access;
  template< typename Archive >
    void serialize( Archive&, const unsigned int /*version*/ );  
};

//
// Traxel datatype
//
class Traxel 
{
 public:
   // construction / assignment
   //takes ownership of locator pointer
  PGMLINK_EXPORT Traxel(unsigned int id = 0,
                        int timestep = 0,
//                        FeatureMap fmap = FeatureMap(),
                        Locator* l = new ComLocator(),
		                ComCorrLocator* lc = new ComCorrLocator()) 
  : Id(id),
    Timestep(timestep),
    locator_(l),
    corr_locator_(lc)
  {
      features = FeatureMapAccessor(this);
  }

  PGMLINK_EXPORT Traxel(const Traxel& other);
  PGMLINK_EXPORT Traxel& operator=(const Traxel& other);
  PGMLINK_EXPORT ~Traxel() { delete locator_; }
  PGMLINK_EXPORT Traxel& set_locator(Locator*);
  PGMLINK_EXPORT Locator* locator() {return locator_;}
  PGMLINK_EXPORT boost::shared_ptr<FeatureStore> get_feature_store() const;
   
   // fields
   unsigned int Id; // id of connected component (aka "label")
   int Timestep; // traxel occured after

   /// before the introduction of TraxelStoreFeatures, the features of a traxel were a member of the traxel,
   /// and one could access the features of a traxel using traxel.features["feat_name"].
   /// FeatureMapAccessor mimics this behaviour in the [] operator, and a traxel gets an instance with the same name.
   class FeatureMapAccessor{
   public:
       FeatureMapAccessor(){}

       /// constructor
       FeatureMapAccessor(Traxel* parent):
           parent_(parent)
       {}

       FeatureMapAccessor(Traxel* parent, const FeatureMap& fm);

       const size_t size() const;

       // non-const access and iterators
       feature_array& operator[](const std::string& feature_name);
       FeatureMap::iterator find(const std::string& feature_name);
       FeatureMap::iterator begin();
       FeatureMap::iterator end();

       // const access and iterators
       FeatureMap::const_iterator find(const std::string& feature_name) const;
       FeatureMap::const_iterator begin() const;
       FeatureMap::const_iterator end() const;
       const FeatureMap& get() const;

       /// This is invoked by a traxel that is added to a feature store
       void feature_store_set_notification();

   private:
       Traxel* parent_;
       FeatureMap feature_map_;

       // boost serialize for Traxel datatype
       friend class boost::serialization::access;
       template< typename Archive >
         void serialize( Archive&, const unsigned int /*version*/ );
   };
   // to allow access via traxel.features["name"]
   FeatureMapAccessor features;

   // method to set the used feature store from outside
   void set_feature_store(boost::shared_ptr<FeatureStore> fs);
   
   // position according to locator
   PGMLINK_EXPORT double X() const;
   PGMLINK_EXPORT double Y() const;
   PGMLINK_EXPORT double Z() const;
    
   PGMLINK_EXPORT double X_corr() const;
   PGMLINK_EXPORT double Y_corr() const;
   PGMLINK_EXPORT double Z_corr() const;

   // relation to other traxels
   PGMLINK_EXPORT double distance_to(const Traxel& other) const;
   PGMLINK_EXPORT double distance_to_corr(const Traxel& other) const;
   PGMLINK_EXPORT double angle(const Traxel& leg1, const Traxel& leg2) const;
   friend std::ostream& operator<< (std::ostream &out, const Traxel &t);


 private:
   boost::shared_ptr<FeatureStore> features_;

   // boost serialize for Traxel datatype
   friend class boost::serialization::access;
   template< typename Archive >
     void serialize( Archive&, const unsigned int /*version*/ );

   Locator* locator_;

   ComCorrLocator* corr_locator_;
 };

 // compare by (time,id) (Traxels can be used as keys (for instance in a std::map) )
 PGMLINK_EXPORT bool operator<(const Traxel& t1, const Traxel& t2);
 PGMLINK_EXPORT bool operator>(const Traxel& t1, const Traxel& t2);
 PGMLINK_EXPORT bool operator==(const Traxel& t1, const Traxel& t2);


 //
 // Traxel collections
 //
 
 // Traxel map
 typedef std::map<unsigned int, Traxel> Traxels;
 template<typename InputIterator>
   Traxels traxel_map_from_traxel_sequence(InputIterator begin, InputIterator end);
 
 // Traxel key-value store
 using namespace boost;
 using namespace boost::multi_index;

 //tags
 struct PGMLINK_EXPORT by_timestep {};
 struct PGMLINK_EXPORT by_timeid {};

 typedef PGMLINK_EXPORT multi_index_container<
   Traxel,
     indexed_by<
       ordered_non_unique<tag<by_timestep>, member<Traxel,int,&Traxel::Timestep> >,
       hashed_unique<
	   tag<by_timeid>,
	   composite_key<
	       Traxel,
	       member<Traxel, int, &Traxel::Timestep>,
	       member<Traxel, unsigned int, &Traxel::Id>
	   > 
	> 
     >
   > 
   TraxelStore;
 typedef PGMLINK_EXPORT TraxelStore::index<by_timestep>::type
   TraxelStoreByTimestep;
 typedef PGMLINK_EXPORT TraxelStore::index<by_timeid>::type
   TraxelStoreByTimeid;

 //
 // TraxelStore functions 
 //
 /**
  * Tight bounding box surrounding the Traxels in store.
  * @return lt lx ly lz ut ux uy uz 
  */
 PGMLINK_EXPORT std::vector<double> bounding_box(const TraxelStore&);

 // timesteps
 PGMLINK_EXPORT std::set<TraxelStoreByTimestep::key_type>
   timesteps(const TraxelStore&);

 PGMLINK_EXPORT TraxelStoreByTimestep::key_type
   earliest_timestep(const TraxelStore&);

 PGMLINK_EXPORT TraxelStoreByTimestep::key_type
   latest_timestep(const TraxelStore&);

 // io
 PGMLINK_EXPORT TraxelStore& add(TraxelStore&, boost::shared_ptr<FeatureStore> fs, Traxel &);
 PGMLINK_EXPORT void set_feature_store(TraxelStore&, boost::shared_ptr<FeatureStore>);

 template<typename InputIt>
   TraxelStore& add(TraxelStore&, InputIt begin, InputIt end);

 std::vector<std::vector<Traxel> > nested_vec_from(const TraxelStore&);

 /** 
  * Filter by field of fiew. 
  * This function adds Traxels from in to out, that are contained in the field of fiew.
  * @return the number of Traxels in the field of view
  */
 class FieldOfView;
 PGMLINK_EXPORT size_t filter_by_fov(const TraxelStore &in, TraxelStore& out, const FieldOfView& );

/**/
/* implementation */
/**/

template< typename Archive >
void Locator::serialize( Archive& ar, const unsigned int /*version*/ ) {
  ar & x_scale;
  ar & y_scale;
  ar & z_scale;
  ar & feature_name_;
}
template< typename Archive >
void ComLocator::serialize( Archive& ar, const unsigned int /*version*/ ) {
  ar & boost::serialization::base_object<Locator>(*this);
}
template< typename Archive >
void IntmaxposLocator::serialize( Archive& ar, const unsigned int /*version*/ ) {
  ar & boost::serialization::base_object<Locator>(*this);
}

template< typename Archive >
void Traxel::serialize( Archive& ar, const unsigned int /*version*/ ) {
  ar.template register_type<ComLocator>();
  ar.template register_type<IntmaxposLocator>();

  ar & Id;
  ar & Timestep;
  ar & features;
  ar & locator_;
}

template< typename Archive >
void Traxel::FeatureMapAccessor::serialize( Archive& ar, const unsigned int /*version*/ ) {
  ar & feature_map_;
  ar & parent_;
}

template<typename InputIterator>
Traxels traxel_map_from_traxel_sequence(InputIterator begin, InputIterator end) {
    Traxels ret;
    for(;begin != end; ++begin) {
	ret[begin->Id] = *begin;
    }
    return ret;
}

template<typename InputIt>
  TraxelStore& add(TraxelStore& ts, InputIt begin, InputIt end) {
  ts.get<by_timestep>().insert(begin, end);
  return ts;
}

} /* namespace pgmlink */


#endif /* TRAXELS_H */
